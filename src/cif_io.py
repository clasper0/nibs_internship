#!/usr/bin/env python3
"""
cif_io.py -- Self-contained CIF parsing, symmetry expansion and structure export.

This module factors out the CIF reader that was previously copy-pasted across
``cif2pdb.py``, ``cif2mols_sdf.py`` and ``cif_to_supercell_sdf.py``.

It implements a minimal subset of the ASE / pymatgen CIF API so that the
conversion pipeline has **no hard dependency on ASE** -- only numpy -- which
makes the COD pre-processing stage easy to deploy on a cluster login node.

Main entry points
-----------------
read_cif(path)                -> Atoms   (full unit cell, symmetry expanded)
make_supercell(atoms, P)      -> Atoms
write_pdb(path, atoms)        -> None
write_sdf(path, atoms)        -> None
cellpar_from_cif(path)        -> (a, b, c, alpha, beta, gamma)

Adapted from the ASE ``cif.py`` reader (LGPL) and trimmed to the tags that
appear in COD / CSD depositions.
"""

from __future__ import annotations


import collections.abc
import re
import shlex
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

CIFDataValue = Union[str, int, float]
CIFData = Union[CIFDataValue, List[CIFDataValue]]


class Atoms:
    """Simple Atoms class to represent atomic structures."""
    
    def __init__(self, symbols=None, positions=None, cell=None, pbc=None):
        self.symbols = symbols or []
        self.positions = np.array(positions) if positions is not None else np.array([])
        self.cell = np.array(cell) if cell is not None else np.zeros((3, 3))
        self.pbc = pbc or [True, True, True]
        self.info = {}
        self.arrays = {}
        
    def __len__(self):
        return len(self.symbols)
    
    def get_chemical_formula(self):
        """Get chemical formula from symbols."""
        if not self.symbols:
            return ""
        
        from collections import Counter
        counts = Counter(self.symbols)
        formula = ""
        for symbol in sorted(counts.keys()):
            count = counts[symbol]
            if count == 1:
                formula += symbol
            else:
                formula += f"{symbol}{count}"
        return formula
    
    def get_scaled_positions(self, wrap=False):
        """Get fractional coordinates."""
        if self.cell.shape == (3, 3) and np.linalg.det(self.cell) != 0:
            inv_cell = np.linalg.inv(self.cell)
            scaled_pos = np.dot(self.positions, inv_cell.T)
            if wrap:
                scaled_pos = scaled_pos % 1.0
            return scaled_pos
        return self.positions.copy()
    
    def set_scaled_positions(self, scaled_positions):
        """Set positions from fractional coordinates."""
        if self.cell.shape == (3, 3) and np.linalg.det(self.cell) != 0:
            self.positions = np.dot(scaled_positions, self.cell)
        else:
            self.positions = np.array(scaled_positions)
    
    def set_array(self, name, array):
        """Set array for atoms."""
        self.arrays[name] = np.array(array)
    
    def wrap(self, eps=1e-7):
        """Wrap atoms into unit cell."""
        scaled = self.get_scaled_positions()
        scaled = scaled % 1.0
        self.set_scaled_positions(scaled)


class Cell:
    """Simple Cell class for unit cell operations."""
    
    def __init__(self, cell_array):
        self.array = np.array(cell_array)
        self.rank = 3 if self.array.shape == (3, 3) and np.linalg.det(self.array) != 0 else 0
    
    def cellpar(self):
        """Get cell parameters [a, b, c, alpha, beta, gamma]."""
        if self.rank != 3:
            return [0, 0, 0, 0, 0, 0]
        
        cell = self.array
        a = np.linalg.norm(cell[0])
        b = np.linalg.norm(cell[1])
        c = np.linalg.norm(cell[2])
        
        alpha = np.arccos(np.dot(cell[1], cell[2]) / (b * c)) * 180 / np.pi
        beta = np.arccos(np.dot(cell[0], cell[2]) / (a * c)) * 180 / np.pi
        gamma = np.arccos(np.dot(cell[0], cell[1]) / (a * b)) * 180 / np.pi
        
        return [a, b, c, alpha, beta, gamma]
    
    @classmethod
    def new(cls, cellpar):
        """Create cell from parameters."""
        if len(cellpar) == 6:
            a, b, c, alpha, beta, gamma = cellpar
            alpha_r = np.radians(alpha)
            beta_r = np.radians(beta)
            gamma_r = np.radians(gamma)
            
            cos_alpha = np.cos(alpha_r)
            cos_beta = np.cos(beta_r)
            cos_gamma = np.cos(gamma_r)
            sin_gamma = np.sin(gamma_r)
            
            cell = np.zeros((3, 3))
            cell[0, 0] = a
            cell[1, 0] = b * cos_gamma
            cell[1, 1] = b * sin_gamma
            cell[2, 0] = c * cos_beta
            cell[2, 1] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
            cell[2, 2] = c * np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 
                                   + 2 * cos_alpha * cos_beta * cos_gamma) / sin_gamma
            return cls(cell)
        else:
            return cls(np.array(cellpar))


def handle_subscripts(text):
    """Handle subscripts in text."""
    return text  # Simplified implementation


def convert_value(value: str) -> CIFDataValue:
    """Convert CIF value string to corresponding python type."""
    value = value.strip()
    if re.match('(".*")|(\'.*\')$', value):
        return handle_subscripts(value[1:-1])
    elif re.match(r'[+-]?\d+$', value):
        return int(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$', value):
        return float(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+\)$', value):
        return float(value[:value.index('(')])
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+$', value):
        warnings.warn(f'Badly formed number: "{value}"')
        return float(value[:value.index('(')])
    else:
        return handle_subscripts(value)


def parse_multiline_string(lines: List[str], line: str) -> str:
    """Parse semicolon-enclosed multiline string and return it."""
    assert line[0] == ';'
    strings = [line[1:].lstrip()]
    while True:
        line = lines.pop().strip()
        if line[:1] == ';':
            break
        strings.append(line)
    return '\n'.join(strings).strip()


def parse_singletag(lines: List[str], line: str) -> Tuple[str, CIFDataValue]:
    """Parse a CIF tag. Returns a key-value pair."""
    kv = line.split(None, 1)
    if len(kv) == 1:
        key = line
        line = lines.pop().strip()
        while not line or line[0] == '#':
            line = lines.pop().strip()
        if line[0] == ';':
            value = parse_multiline_string(lines, line)
        else:
            value = line
    else:
        key, value = kv
    return key, convert_value(value)


def parse_cif_loop_headers(lines: List[str]) -> Iterator[str]:
    while lines:
        line = lines.pop()
        tokens = line.split()

        if len(tokens) == 1 and tokens[0].startswith('_'):
            header = tokens[0].lower()
            yield header
        else:
            lines.append(line)
            return


def parse_cif_loop_data(lines: List[str], ncolumns: int) -> List[List[CIFDataValue]]:
    columns: List[List[CIFDataValue]] = [[] for _ in range(ncolumns)]

    tokens = []
    while lines:
        line = lines.pop().strip()
        lowerline = line.lower()
        if (not line or
            line.startswith('_') or
            lowerline.startswith('data_') or
                lowerline.startswith('loop_')):
            lines.append(line)
            break

        if line.startswith('#'):
            continue

        line = line.split(' #')[0]

        if line.startswith(';'):
            moretokens = [parse_multiline_string(lines, line)]
        else:
            if ncolumns == 1:
                moretokens = [line]
            else:
                moretokens = shlex.split(line, posix=False)

        tokens += moretokens
        if len(tokens) < ncolumns:
            continue
        if len(tokens) == ncolumns:
            for i, token in enumerate(tokens):
                columns[i].append(convert_value(token))
        else:
            warnings.warn(f'Wrong number {len(tokens)} of tokens, '
                          f'expected {ncolumns}: {tokens}')

        tokens = []

    return columns


def parse_loop(lines: List[str]) -> Dict[str, List[CIFDataValue]]:
    """Parse a CIF loop. Returns a dict with column tag names as keys."""
    headers = list(parse_cif_loop_headers(lines))
    columns = parse_cif_loop_data(lines, len(headers))

    columns_dict = {}
    for i, header in enumerate(headers):
        if header in columns_dict:
            warnings.warn(f'Duplicated loop tags: {header}')
        else:
            columns_dict[header] = columns[i]
    return columns_dict


def parse_items(lines: List[str], line: str) -> Dict[str, CIFData]:
    """Parse CIF data items and return a dict with all tags."""
    tags: Dict[str, CIFData] = {}

    while True:
        if not lines:
            break
        line = lines.pop().strip()
        if not line:
            continue
        lowerline = line.lower()
        if not line or line.startswith('#'):
            continue
        elif line.startswith('_'):
            key, value = parse_singletag(lines, line)
            tags[key.lower()] = value
        elif lowerline.startswith('loop_'):
            tags.update(parse_loop(lines))
        elif lowerline.startswith('data_'):
            if line:
                lines.append(line)
            break
        elif line.startswith(';'):
            parse_multiline_string(lines, line)
        else:
            raise ValueError(f'Unexpected CIF file entry: "{line}"')
    return tags


def parse_symmetry_operation(symop_str):
    """Parse a symmetry operation string like 'x,y,z' or 'x+1/2,-y,z'."""
    # Remove spaces and split by commas
    parts = [p.strip() for p in symop_str.replace(' ', '').split(',')]
    
    if len(parts) != 3:
        raise ValueError(f"Invalid symmetry operation: {symop_str}")
    
    # Initialize rotation matrix and translation vector
    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    
    for i, part in enumerate(parts):
        # Handle each coordinate (x, y, z)
        for j, var in enumerate(['x', 'y', 'z']):
            if var in part:
                # Find coefficient
                var_idx = part.find(var)
                if var_idx == 0 or part[var_idx-1] == '+':
                    coeff = 1.0
                elif part[var_idx-1] == '-':
                    coeff = -1.0
                else:
                    # Extract coefficient
                    start = var_idx - 1
                    while start > 0 and part[start-1] not in '+-':
                        start -= 1
                    coeff_str = part[start:var_idx]
                    if coeff_str == '' or coeff_str == '+':
                        coeff = 1.0
                    elif coeff_str == '-':
                        coeff = -1.0
                    else:
                        coeff = float(coeff_str)
                
                rotation[i, j] = coeff
        
        # Handle translation part
        # Remove variables to get translation terms
        trans_part = part
        for var in ['x', 'y', 'z']:
            trans_part = re.sub(r'[+-]?[\d]*' + var, '', trans_part)
        
        # Parse remaining translation
        if trans_part:
            # Handle fractions and decimals
            trans_value = 0.0
            # Split by + and - while keeping the signs
            terms = re.split(r'(?=[+-])', trans_part)
            for term in terms:
                if term.strip():
                    if '/' in term:
                        # Handle fractions
                        if '+' in term or '-' in term:
                            sign = 1 if term.startswith('+') else -1 if term.startswith('-') else 1
                            frac_part = term.lstrip('+-')
                        else:
                            sign = 1
                            frac_part = term
                        
                        if '/' in frac_part:
                            num, den = frac_part.split('/')
                            trans_value += sign * float(num) / float(den)
                        else:
                            trans_value += sign * float(frac_part)
                    else:
                        # Handle decimals
                        if term.strip():
                            trans_value += float(term)
            
            translation[i] = trans_value
    
    return rotation, translation


class CIFBlock:
    """A block in a crystallographic information file."""

    cell_tags = ['_cell_length_a', '_cell_length_b', '_cell_length_c',
                 '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']

    def __init__(self, name: str, tags: Dict[str, CIFData]):
        self.name = name
        self._tags = tags

    def __getitem__(self, key: str) -> CIFData:
        return self._tags[key]

    def __contains__(self, key):
        return key in self._tags

    def get(self, key, default=None):
        return self._tags.get(key, default)

    def get_cellpar(self) -> Optional[List]:
        try:
            return [self[tag] for tag in self.cell_tags]
        except KeyError:
            return None

    def get_cell(self) -> Cell:
        cellpar = self.get_cellpar()
        if cellpar is None:
            return Cell(np.zeros((3, 3)))
        return Cell.new(cellpar)

    def _raw_scaled_positions(self) -> Optional[np.ndarray]:
        coords = [self.get(name) for name in ['_atom_site_fract_x',
                                              '_atom_site_fract_y',
                                              '_atom_site_fract_z']]
        if None in coords:
            return None
        return np.array(coords).T

    def get_symbols(self) -> List[str]:
        labels = self._get_any(['_atom_site_type_symbol', '_atom_site_label'])
        if labels is None:
            raise ValueError('No symbols found')

        symbols = []
        for label in labels:
            if label == '.' or label == '?':
                raise ValueError('Symbols are undetermined')
            # Strip off additional labeling on chemical symbols
            match = re.search(r'([A-Z][a-z]?)', str(label))
            if match:
                symbol = match.group(0)
                # Replace deuterium with hydrogen
                if symbol == 'D':
                    symbol = 'H'
                symbols.append(symbol)
            else:
                symbols.append('X')  # Unknown element
        return symbols

    def _get_any(self, names):
        for name in names:
            if name in self:
                return self[name]
        return None

    def get_symmetry_operations(self):
        """Get symmetry operations from CIF."""
        # Try different possible symmetry operation tags
        symop_tags = [
            '_symmetry_equiv_pos_as_xyz',
            '_space_group_symop_operation_xyz',
            '_symmetry_equiv_pos_xyz'
        ]
        
        symops = None
        for tag in symop_tags:
            symops = self.get(tag)
            if symops is not None:
                break
        
        if symops is None:
            # Default to identity operation if no symmetry operations found
            print("Warning: No symmetry operations found, using identity operation")
            return [np.eye(3), np.zeros(3)]
        
        if isinstance(symops, str):
            symops = [symops]
        
        operations = []
        for symop in symops:
            try:
                rotation, translation = parse_symmetry_operation(symop)
                operations.append((rotation, translation))
            except Exception as e:
                print(f"Warning: Could not parse symmetry operation '{symop}': {e}")
                continue
        
        if not operations:
            # Fallback to identity
            operations = [(np.eye(3), np.zeros(3))]
        
        return operations

    def get_occupancies(self):
        """Get site occupancies."""
        occs = self.get('_atom_site_occupancy')
        if occs is None:
            # If no occupancy data, assume all sites are fully occupied
            symbols = self.get_symbols()
            return [1.0] * len(symbols)
        return occs

    def get_atoms(self) -> Atoms:
        """Returns an Atoms object from CIF tags with all symmetry-generated atoms."""
        cell = self.get_cell()
        
        # Get atomic symbols
        base_symbols = self.get_symbols()
        
        # Get positions
        base_positions = self._raw_scaled_positions()
        if base_positions is None:
            raise ValueError('No atomic positions found')
        
        # Get occupancies
        occupancies = self.get_occupancies()
        
        # Get symmetry operations
        symops = self.get_symmetry_operations()
        
        print(f"Found {len(base_symbols)} unique atomic sites")
        print(f"Found {len(symops)} symmetry operations")
        
        # Apply symmetry operations to generate all atoms in the unit cell
        all_symbols = []
        all_positions = []
        tolerance = 1e-5
        
        for rotation, translation in symops:
            for i, (symbol, pos, occ) in enumerate(zip(base_symbols, base_positions, occupancies)):
                # Skip sites with zero occupancy
                if abs(occ) < tolerance:
                    continue
                
                # Apply symmetry operation
                new_pos = np.dot(rotation, pos) + translation
                
                # Wrap to unit cell
                new_pos = new_pos % 1.0
                
                # Check if this position already exists (avoid duplicates)
                is_duplicate = False
                for existing_pos in all_positions:
                    diff = new_pos - existing_pos
                    # Handle periodic boundary conditions
                    diff = diff - np.round(diff)
                    if np.linalg.norm(diff) < tolerance:
                        is_duplicate = True
                        break
                
                if not is_duplicate:
                    all_symbols.append(symbol)
                    all_positions.append(new_pos)
        
        print(f"Generated {len(all_symbols)} atoms in the complete unit cell")
        
        # Create atoms object
        atoms = Atoms(symbols=all_symbols, cell=cell.array)
        atoms.set_scaled_positions(np.array(all_positions))
        
        return atoms


def parse_block(lines: List[str], line: str) -> CIFBlock:
    assert line.lower().startswith('data_')
    blockname = line.split('_', 1)[1].rstrip()
    tags = parse_items(lines, line)
    return CIFBlock(blockname, tags)


def parse_cif(fileobj) -> Iterator[CIFBlock]:
    """Parse a CIF file."""
    if isinstance(fileobj, (str, Path)):
        with open(fileobj, 'rb') as f:
            data = f.read()
    else:
        data = fileobj.read()

    if isinstance(data, bytes):
        data = data.decode('latin1')
    
    lines = [e for e in data.split('\n') if len(e) > 0]
    lines = [''] + lines[::-1]

    while lines:
        line = lines.pop().strip()
        if not line or line.startswith('#'):
            continue
        yield parse_block(lines, line)


def read_cif(filename):
    """Read CIF file and return Atoms object."""
    blocks = list(parse_cif(filename))
    if not blocks:
        raise ValueError("No data blocks found in CIF file")
    
    # Use first block with structure data
    for block in blocks:
        try:
            atoms = block.get_atoms()
            return atoms
        except (ValueError, KeyError) as e:
            print(f"Warning: Skipping block {block.name}: {e}")
            continue
    
    raise ValueError("No valid structure data found in CIF file")


def make_supercell(prim, P, wrap=True, order="cell-major", tol=1e-5):
    """Generate a supercell by applying transformation matrix P."""
    
    def lattice_points_in_supercell(supercell_matrix):
        """Find all lattice points contained in a supercell."""
        diagonals = np.array([
            [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
            [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1],
        ])
        d_points = np.dot(diagonals, supercell_matrix)

        mins = np.min(d_points, axis=0)
        maxes = np.max(d_points, axis=0) + 1

        ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
        br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
        cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

        all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
        all_points = all_points.reshape((-1, 3))

        frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))
        tvects = frac_points[np.all(frac_points < 1 - 1e-10, axis=1)
                             & np.all(frac_points >= -1e-10, axis=1)]
        return tvects

    supercell_matrix = P
    supercell_cell = supercell_matrix @ prim.cell

    # Get lattice points
    lattice_points_frac = lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell_cell)
    N = len(lattice_points)

    if order == "cell-major":
        shifted = prim.positions[None, :, :] + lattice_points[:, None, :]
    elif order == "atom-major":
        shifted = prim.positions[:, None, :] + lattice_points[None, :, :]
    else:
        raise ValueError(f"invalid order: {order}")
    
    shifted_reshaped = shifted.reshape(-1, 3)

    # Create supercell atoms
    superatoms = Atoms(
        symbols=prim.symbols * N,
        positions=shifted_reshaped,
        cell=supercell_cell,
        pbc=prim.pbc
    )

    # Copy arrays
    for name, arr in prim.arrays.items():
        if name == "positions":
            continue
        shape = (N * arr.shape[0], *arr.shape[1:])
        if order == "cell-major":
            new_arr = np.repeat(arr[None, :], N, axis=0).reshape(shape)
        elif order == "atom-major":
            new_arr = np.repeat(arr[:, None], N, axis=1).reshape(shape)
        superatoms.set_array(name, new_arr)

    if wrap:
        superatoms.wrap(eps=tol)

    return superatoms


def write_pdb(filename, atoms):
    """Write atoms to PDB format."""
    
    # Element to radius mapping for PDB (rough approximation)
    atomic_radii = {
        'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
        'P': 1.8, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
    }
    
    with open(filename, 'w') as f:
        f.write("TITLE     Generated supercell structure from complete unit cell\n")
        f.write("REMARK    Created by CIF to Supercell PDB converter\n")
        f.write("REMARK    Includes all symmetry-generated atoms\n")
        
        # Write CRYST1 record with cell parameters
        if atoms.cell is not None and atoms.cell.shape == (3, 3):
            cell_obj = Cell(atoms.cell)
            cellpar = cell_obj.cellpar()
            f.write(f"CRYST1{cellpar[0]:9.3f}{cellpar[1]:9.3f}{cellpar[2]:9.3f}"
                   f"{cellpar[3]:7.2f}{cellpar[4]:7.2f}{cellpar[5]:7.2f} P 1           1\n")
        
        # Write atoms
        for i, (symbol, pos) in enumerate(zip(atoms.symbols, atoms.positions)):
            atom_name = f"{symbol:<4}"
            res_name = "UNL"  # Unknown ligand
            chain_id = "A"
            res_seq = 1
            
            f.write(f"HETATM{i+1:5d} {atom_name} {res_name} {chain_id}{res_seq:4d}    "
                   f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                   f"{1.00:6.2f}{0.00:6.2f}          {symbol:>2}\n")
        
        f.write("END\n")




def write_sdf(filename, atoms, title=None):
    """Write atoms to a V2000 SDF file.

    Only coordinates and element symbols are emitted -- bond perception is
    deliberately left to Unicon / RDKit downstream, which is what the COD
    pipeline relies on.
    """
    title = title or "CIF-derived structure"
    symbols = list(atoms.symbols)
    positions = np.asarray(atoms.positions)
    n = len(symbols)

    with open(filename, "w") as f:
        f.write(f"{title}\n")
        f.write("  cif_io  COD/CSD pre-processing\n")
        f.write("\n")
        f.write(f"{n:3d}{0:3d}  0  0  0  0  0  0  0  0999 V2000\n")
        for sym, pos in zip(symbols, positions):
            f.write(f"{pos[0]:10.4f}{pos[1]:10.4f}{pos[2]:10.4f} "
                    f"{sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n")
        f.write("M  END\n")
        f.write("$$$$\n")


def cellpar_from_cif(path):
    """Return (a, b, c, alpha, beta, gamma) read straight from the CIF header.

    Returns ``None`` when any of the six cell tags is absent.
    """
    tags = {
        "_cell_length_a": None,
        "_cell_length_b": None,
        "_cell_length_c": None,
        "_cell_angle_alpha": None,
        "_cell_angle_beta": None,
        "_cell_angle_gamma": None,
    }
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.strip()
            for tag in tags:
                if stripped.startswith(tag):
                    value = stripped[len(tag):].strip()
                    value = value.split("(")[0]          # drop e.s.d. e.g. 7.8783(2)
                    try:
                        tags[tag] = float(value)
                    except ValueError:
                        pass
    if any(v is None for v in tags.values()):
        return None
    return tuple(tags[t] for t in
                 ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                  "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))


def supercell_from_cif(path, n=1):
    """Convenience: read a CIF and expand it n x n x n."""
    atoms = read_cif(path)
    if n != 1:
        atoms = make_supercell(atoms, np.array([[n, 0, 0], [0, n, 0], [0, 0, n]]))
    return atoms


def _main(argv=None):
    import argparse

    ap = argparse.ArgumentParser(description="CIF -> PDB / SDF structure converter")
    ap.add_argument("input_cif", help="input .cif file")
    ap.add_argument("-o", "--output", help="output .pdb or .sdf (default: alongside input)")
    ap.add_argument("-n", "--supercell", type=int, default=1,
                    help="expand to n x n x n supercell (default: 1 = unit cell only)")
    ap.add_argument("--cellpar", action="store_true",
                    help="only print the cell parameters and exit")
    args = ap.parse_args(argv)

    if args.cellpar:
        print(cellpar_from_cif(args.input_cif))
        return 0

    atoms = supercell_from_cif(args.input_cif, args.supercell)
    out = args.output or str(Path(args.input_cif).with_suffix(".pdb"))

    print(f"Loaded {args.input_cif}")
    print(f"  atoms : {len(atoms)}")
    print(f"  formula: {atoms.get_chemical_formula()}")

    if out.lower().endswith(".sdf"):
        write_sdf(out, atoms)
    else:
        write_pdb(out, atoms)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
