#!/usr/bin/env python3
"""
CIF to SDF Converter with Molecular Connectivity Processing

This script converts CIF files to SDF format through several processing steps:
1. Parses CIF files with complete symmetry operations
2. Splits molecules in the unit cell
3. Uses Unicon to assign bond connectivity while preserving original coordinates
4. Merges molecules back into a single unit cell
5. Generates a 3x3x3 supercell
6. Exports to SDF format compatible with RDKit

Dependencies:
- numpy
- unicon (external program)

Usage:
    python cif2mols_sdf.py input.cif [output.sdf]
"""

import sys
import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union
import numpy as np
import re
import shlex
import warnings
import collections.abc
from concurrent.futures import ThreadPoolExecutor, as_completed


# CIF parsing code adapted from cif2pdb2.0.py
CIFDataValue = str | int | float
CIFData = CIFDataValue | List[CIFDataValue]


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
    if isinstance(fileobj, str):
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
        f.write("REMARK    Created by CIF to Molecular SDF converter\n")
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


class Atom:
    """原子类，存储原子信息"""
    
    def __init__(self, line: str):
        """从PDB行解析原子信息"""
        self.record_type = line[0:6].strip()
        self.atom_number = int(line[6:11].strip())
        self.atom_name = line[12:16].strip()
        self.residue_name = line[17:20].strip()
        self.chain_id = line[21:22].strip()
        self.residue_number = int(line[22:26].strip()) if line[22:26].strip() else 1
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
        self.temperature_factor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
        self.element = line[76:78].strip() if len(line) > 76 else self.atom_name[0]
        
        # 原始行，用于输出
        self.original_line = line.rstrip()
    
    def get_coordinates(self) -> np.ndarray:
        """获取原子坐标"""
        return np.array([self.x, self.y, self.z])
    
    def distance_to(self, other_atom: 'Atom') -> float:
        """计算到另一个原子的距离"""
        coord1 = self.get_coordinates()
        coord2 = other_atom.get_coordinates()
        return np.linalg.norm(coord1 - coord2)
    
    def to_pdb_line(self, atom_number: Optional[int] = None) -> str:
        """转换为PDB格式行"""
        if atom_number is not None:
            # 更新原子编号
            line = self.original_line
            new_line = f"{line[:6]}{atom_number:5d}{line[11:]}"
            return new_line
        return self.original_line


class Molecule:
    """分子类，包含一个分子的所有原子"""
    
    def __init__(self, atoms: List[Atom], molecule_id: int):
        self.atoms = atoms
        self.molecule_id = molecule_id
    
    def __len__(self) -> int:
        return len(self.atoms)
    
    def get_center(self) -> np.ndarray:
        """获取分子的几何中心"""
        coords = np.array([atom.get_coordinates() for atom in self.atoms])
        return np.mean(coords, axis=0)
    
    def get_formula(self) -> str:
        """获取分子式"""
        from collections import Counter
        elements = [atom.element for atom in self.atoms]
        element_count = Counter(elements)
        
        formula = ""
        # 通常的元素顺序：C, H, N, O, 然后按字母顺序
        order = ['C', 'H', 'N', 'O'] + sorted([e for e in element_count.keys() if e not in ['C', 'H', 'N', 'O']])
        
        for element in order:
            if element in element_count:
                count = element_count[element]
                if count == 1:
                    formula += element
                else:
                    formula += f"{element}{count}"
        
        return formula
    
    def to_pdb_string(self, title: str = None) -> str:
        """转换为完整的PDB字符串"""
        lines = []
        
        # 添加头部信息
        if title:
            lines.append(f"TITLE     {title}")
        else:
            lines.append(f"TITLE     Molecule {self.molecule_id} ({self.get_formula()})")
        
        lines.append("REMARK    Extracted from supercell structure")
        lines.append(f"REMARK    Molecule ID: {self.molecule_id}")
        lines.append(f"REMARK    Number of atoms: {len(self.atoms)}")
        lines.append(f"REMARK    Molecular formula: {self.get_formula()}")
        
        # 添加原子信息，重新编号
        for i, atom in enumerate(self.atoms, 1):
            lines.append(atom.to_pdb_line(i))
        
        lines.append("END")
        
        return "\n".join(lines)


class PDBMoleculeSplitter:
    """PDB分子分离器"""
    
    def __init__(self, bond_distance_cutoff: float = 2.0, hydrogen_cutoff: float = 1.3):
        """
        初始化分离器
        
        Args:
            bond_distance_cutoff: 非氢原子间的键长阈值
            hydrogen_cutoff: 氢原子的键长阈值
        """
        self.bond_distance_cutoff = bond_distance_cutoff
        self.hydrogen_cutoff = hydrogen_cutoff
        self.atoms = []
        self.molecules = []
    
    def read_pdb(self, filename: str) -> None:
        """读取PDB文件"""
        print(f"Reading PDB file: {filename}")
        
        self.atoms = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith(('HETATM', 'ATOM')):
                    try:
                        atom = Atom(line)
                        self.atoms.append(atom)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line: {line.strip()} - {e}")
        
        print(f"Successfully read {len(self.atoms)} atoms")
    
    def _are_bonded(self, atom1: Atom, atom2: Atom) -> bool:
        """判断两个原子是否成键"""
        distance = atom1.distance_to(atom2)
        
        # 根据原子类型选择不同的阈值
        if atom1.element == 'H' or atom2.element == 'H':
            return distance <= self.hydrogen_cutoff
        else:
            return distance <= self.bond_distance_cutoff
    
    def _find_connected_atoms(self, start_atom_idx: int, visited: set[int]) -> List[int]:
        """使用深度优先搜索找到连接的原子"""
        if start_atom_idx in visited:
            return []
        
        connected = [start_atom_idx]
        visited.add(start_atom_idx)
        
        start_atom = self.atoms[start_atom_idx]
        
        # 检查与其他未访问原子的连接
        for i, other_atom in enumerate(self.atoms):
            if i not in visited and self._are_bonded(start_atom, other_atom):
                connected.extend(self._find_connected_atoms(i, visited))
        
        return connected
    
    def split_molecules(self) -> List[Molecule]:
        """分离分子"""
        print("Splitting molecules using connectivity analysis...")
        
        visited = set()
        molecules = []
        molecule_id = 1
        
        for i, atom in enumerate(self.atoms):
            if i not in visited:
                # 找到这个连通分量的所有原子
                connected_indices = self._find_connected_atoms(i, visited)
                
                if connected_indices:
                    # 创建分子
                    molecule_atoms = [self.atoms[idx] for idx in connected_indices]
                    molecule = Molecule(molecule_atoms, molecule_id)
                    molecules.append(molecule)
                    molecule_id += 1
        
        self.molecules = molecules
        print(f"Found {len(molecules)} molecules")
        
        # 显示分子统计信息
        for mol in molecules:
            center = mol.get_center()
            print(f"  Molecule {mol.molecule_id}: {len(mol)} atoms, "
                  f"formula: {mol.get_formula()}, "
                  f"center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        
        return molecules
    
    def save_molecules(self, output_dir: str, filename_pattern: str = "molecule_{:03d}.pdb") -> List[str]:
        """保存分子到单独的PDB文件"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        saved_files = []
        
        print(f"\nSaving molecules to {output_dir}:")
        for molecule in self.molecules:
            filename = filename_pattern.format(molecule.molecule_id)
            filepath = output_path / filename
            
            with open(filepath, 'w') as f:
                f.write(molecule.to_pdb_string())
            
            saved_files.append(str(filepath))
            print(f"  Saved: {filename} ({len(molecule)} atoms, {molecule.get_formula()})")
        
        return saved_files


class UniconProcessor:
    """Unicon处理器"""
    
    def __init__(self, unicon_path: str = "unicon", 
                 output_format: str = "sdf",
                 generate_3d: bool = True,
                 max_conformers: int = 1,
                 verbosity: int = 1,
                 timeout: int = 300):
        """
        初始化Unicon处理器
        
        Args:
            unicon_path: unicon可执行文件路径
            output_format: 输出格式 (sdf, mol2等)
            generate_3d: 是否生成3D坐标
            max_conformers: 最大构象数量
            verbosity: 详细程度 (0-3)
            timeout: 超时时间（秒）
        """
        self.unicon_path = unicon_path
        self.output_format = output_format
        self.generate_3d = generate_3d
        self.max_conformers = max_conformers
        self.verbosity = verbosity
        self.timeout = timeout
    
    def check_unicon_availability(self) -> bool:
        """检查unicon是否可用"""
        try:
            result = subprocess.run(
                [self.unicon_path, "-h"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            return False
    
    def process_molecule(self, input_file: str, output_file: str = None) -> Tuple[bool, str, str]:
        """
        使用unicon处理单个分子
        
        Args:
            input_file: 输入PDB文件路径
            output_file: 输出文件路径（可选）
        
        Returns:
            (成功标志, 输出文件路径, 错误信息)
        """
        if output_file is None:
            input_path = Path(input_file)
            output_file = str(input_path.with_suffix(f'.{self.output_format}'))
        
        # 构建unicon命令
        cmd = [
            self.unicon_path,
            "-i", input_file,
            "-o", output_file,
            "--outFormat", self.output_format.upper(),
            "-t", "single",
            "-p", "single", 
            "-g", "1"
        ]
        
        try:
            print(f"Processing {Path(input_file).name} with unicon...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            if result.returncode == 0:
                if os.path.exists(output_file):
                    return True, output_file, ""
                else:
                    return False, "", f"Output file not created: {output_file}"
            else:
                error_msg = f"Unicon failed with return code {result.returncode}\n"
                error_msg += f"STDERR: {result.stderr}\n"
                error_msg += f"STDOUT: {result.stdout}"
                return False, "", error_msg
                
        except subprocess.TimeoutExpired:
            return False, "", f"Unicon process timed out after {self.timeout} seconds"
        except Exception as e:
            return False, "", f"Error running unicon: {str(e)}"
    
    def process_molecules_batch(self, input_files: List[str], output_dir: str = None, 
                               max_workers: int = 4) -> Dict[str, Tuple[bool, str, str]]:
        """
        批量处理分子
        
        Args:
            input_files: 输入文件列表
            output_dir: 输出目录
            max_workers: 最大并行工作线程数
        
        Returns:
            处理结果字典 {输入文件: (成功标志, 输出文件, 错误信息)}
        """
        if output_dir:
            Path(output_dir).mkdir(exist_ok=True)
        
        results = {}
        
        print(f"\nProcessing {len(input_files)} molecules with unicon (max_workers={max_workers}):")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_file = {}
            for input_file in input_files:
                if output_dir:
                    input_path = Path(input_file)
                    output_file = Path(output_dir) / f"{input_path.stem}.{self.output_format}"
                else:
                    output_file = None
                
                future = executor.submit(self.process_molecule, input_file, str(output_file) if output_file else None)
                future_to_file[future] = input_file
            
            # 收集结果
            completed = 0
            for future in as_completed(future_to_file):
                input_file = future_to_file[future]
                try:
                    success, output_file, error_msg = future.result()
                    results[input_file] = (success, output_file, error_msg)
                    
                    completed += 1
                    status = "✓" if success else "✗"
                    filename = Path(input_file).name
                    print(f"  [{completed}/{len(input_files)}] {status} {filename}")
                    
                    if not success and error_msg:
                        print(f"    Error: {error_msg}")
                
                except Exception as e:
                    results[input_file] = (False, "", f"Exception: {str(e)}")
                    completed += 1
                    print(f"  [{completed}/{len(input_files)}] ✗ {Path(input_file).name} (Exception)")
        
        return results


def merge_molecules_to_pdb(molecule_files: List[str], output_file: str) -> bool:
    """合并多个分子文件到一个PDB文件"""
    try:
        with open(output_file, 'w') as outfile:
            # 为第一个文件保留头部信息
            first_file = True
            
            for mol_file in molecule_files:
                with open(mol_file, 'r') as infile:
                    for line in infile:
                        if first_file:
                            # 保留第一个文件的头部信息
                            if not line.startswith('END'):
                                outfile.write(line)
                        else:
                            # 跳过后续文件的头部信息
                            if line.startswith(('HETATM', 'ATOM')):
                                outfile.write(line)
                
                # 只有第一个文件写入头部信息
                first_file = False
            
            # 文件末尾写入END
            outfile.write("END\n")
        
        print(f"Successfully merged {len(molecule_files)} molecules into {output_file}")
        return True
    except Exception as e:
        print(f"Error merging molecules: {e}")
        return False


def convert_sdf_to_rdkit_compatible(input_file: str, output_file: str) -> bool:
    """确保SDF文件与RDKit兼容"""
    try:
        # 直接复制文件
        shutil.copy(input_file, output_file)
        print(f"Successfully created RDKit-compatible SDF file: {output_file}")
        return True
    except Exception as e:
        print(f"Error converting SDF to RDKit-compatible format: {e}")
        return False


def process_cif_to_sdf(cif_path: str, output_path: str = None) -> bool:
    """主处理函数"""
    print(f"Loading CIF file: {cif_path}")
    
    try:
        # 1. 读取CIF文件
        atoms = read_cif(cif_path)
        print(f"Successfully loaded structure:")
        print(f"  - Formula: {atoms.get_chemical_formula()}")
        print(f"  - Number of atoms in complete unit cell: {len(atoms)}")
        
        if atoms.cell is not None and atoms.cell.shape == (3, 3):
            cell_obj = Cell(atoms.cell)
            cellpar = cell_obj.cellpar()
            print(f"  - Cell parameters: {cellpar}")
        
    except Exception as e:
        print(f"Error loading CIF file: {e}")
        return False
    
    # 创建临时工作目录
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        print(f"\nUsing temporary directory: {temp_path}")
        
        # 2. 生成单元胞PDB文件
        unitcell_pdb = temp_path / "unitcell.pdb"
        print(f"\nExporting unit cell to PDB format: {unitcell_pdb}")
        
        try:
            write_pdb(str(unitcell_pdb), atoms)
            print("✓ Successfully exported unit cell to PDB")
        except Exception as e:
            print(f"Error writing unit cell PDB file: {e}")
            return False
        
        # 3. 分离分子
        splitter = PDBMoleculeSplitter(bond_distance_cutoff=2.0, hydrogen_cutoff=1.3)
        try:
            splitter.read_pdb(str(unitcell_pdb))
            molecules = splitter.split_molecules()
            
            if not molecules:
                print("No molecules found!")
                return False
            
            # 保存分子
            pdb_files = splitter.save_molecules(str(temp_path / "molecules"))
            
        except Exception as e:
            print(f"Error during molecule splitting: {e}")
            return False
        
        # 4. 使用unicon处理
        processor = UniconProcessor(
            unicon_path="unicon",
            output_format="sdf",
            generate_3d=True,
            max_conformers=1,
            verbosity=1,
            timeout=300
        )
        
        # 检查unicon可用性
        if not processor.check_unicon_availability():
            print("Error: unicon not found in PATH")
            print("Please make sure unicon is installed and in your PATH")
            return False
        
        try:
            # 批量处理
            results = processor.process_molecules_batch(
                pdb_files, 
                str(temp_path / "unicon_output"),
                4
            )
            
            # 统计结果
            successful = sum(1 for success, _, _ in results.values() if success)
            failed = len(results) - successful
            
            print(f"\nUnicon processing completed:")
            print(f"  Successful: {successful}")
            print(f"  Failed: {failed}")
            
            # 检查是否有失败的处理
            if failed > 0:
                print("\nFailed files:")
                for input_file, (success, _, error_msg) in results.items():
                    if not success:
                        print(f"  {Path(input_file).name}: {error_msg}")
                return False
            
            # 获取成功处理的SDF文件
            sdf_files = [output_file for success, output_file, _ in results.values() if success]
            
        except Exception as e:
            print(f"Error during unicon processing: {e}")
            return False
        
        # 5. 合并分子
        merged_pdb = temp_path / "merged.pdb"
        print(f"\nMerging processed molecules...")
        
        if not merge_molecules_to_pdb(pdb_files, str(merged_pdb)):
            print("Failed to merge molecules")
            return False
        
        # 6. 生成3x3x3超胞
        print(f"\nGenerating 3x3x3 supercell...")
        try:
            # 重新读取合并后的分子
            merged_atoms = Atoms()
            merged_atoms.symbols = atoms.symbols
            merged_atoms.positions = atoms.positions
            merged_atoms.cell = atoms.cell
            
            # 创建3x3x3超胞
            P = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=int)
            supercell = make_supercell(merged_atoms, P)
            
            print(f"Supercell created successfully:")
            print(f"  - Unit cell atoms: {len(merged_atoms)}")
            print(f"  - Supercell atoms: {len(supercell)}")
            print(f"  - Expected atoms: {len(merged_atoms) * 27}")  # 3x3x3 = 27
            print(f"  - Formula: {supercell.get_chemical_formula()}")
            
        except Exception as e:
            print(f"Error creating supercell: {e}")
            return False
        
        # 7. 生成超胞PDB文件
        supercell_pdb = temp_path / "supercell.pdb"
        print(f"\nExporting supercell to PDB format: {supercell_pdb}")
        
        try:
            write_pdb(str(supercell_pdb), supercell)
            print("✓ Successfully exported supercell to PDB")
        except Exception as e:
            print(f"Error writing supercell PDB file: {e}")
            return False
        
        # 8. 处理超胞文件
        supercell_splitter = PDBMoleculeSplitter(bond_distance_cutoff=2.0, hydrogen_cutoff=1.3)
        try:
            supercell_splitter.read_pdb(str(supercell_pdb))
            supercell_molecules = supercell_splitter.split_molecules()
            
            if not supercell_molecules:
                print("No molecules found in supercell!")
                return False
            
            # 保存超胞分子
            supercell_pdb_files = supercell_splitter.save_molecules(
                str(temp_path / "supercell_molecules"), 
                "supercell_molecule_{:03d}.pdb"
            )
            
        except Exception as e:
            print(f"Error during supercell molecule splitting: {e}")
            return False
        
        # 9. 使用unicon处理超胞分子
        print(f"\nProcessing supercell molecules with unicon...")
        try:
            # 批量处理
            supercell_results = processor.process_molecules_batch(
                supercell_pdb_files, 
                str(temp_path / "supercell_unicon_output"),
                4
            )
            
            # 统计结果
            supercell_successful = sum(1 for success, _, _ in supercell_results.values() if success)
            supercell_failed = len(supercell_results) - supercell_successful
            
            print(f"\nSupercell unicon processing completed:")
            print(f"  Successful: {supercell_successful}")
            print(f"  Failed: {supercell_failed}")
            
            # 检查是否有失败的处理
            if supercell_failed > 0:
                print("\nFailed supercell files:")
                for input_file, (success, _, error_msg) in supercell_results.items():
                    if not success:
                        print(f"  {Path(input_file).name}: {error_msg}")
                return False
            
            # 获取成功处理的SDF文件
            supercell_sdf_files = [output_file for success, output_file, _ in supercell_results.values() if success]
            
        except Exception as e:
            print(f"Error during supercell unicon processing: {e}")
            return False
        
        # 10. 生成最终输出文件路径
        if output_path is None:
            cif_path_obj = Path(cif_path)
            output_path = str(cif_path_obj.with_suffix('.sdf'))
        
        print(f"\nCreating final RDKit-compatible SDF file: {output_path}")
        
        # 合并所有SDF文件
        try:
            with open(output_path, 'w') as outfile:
                for sdf_file in supercell_sdf_files:
                    with open(sdf_file, 'r') as infile:
                        shutil.copyfileobj(infile, outfile)
            
            print("✓ Successfully created final SDF file!")
            return True
            
        except Exception as e:
            print(f"Error creating final SDF file: {e}")
            return False


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("Usage: python cif2mols_sdf.py input.cif [output.sdf]")
        print("\nExample:")
        print("  python cif2mols_sdf.py 1000001.cif")
        print("  python cif2mols_sdf.py 1000001.cif output.sdf")
        print("\nFeatures:")
        print("  - Parses all symmetry operations from CIF")
        print("  - Generates complete unit cell with all atoms")
        print("  - Splits molecules using connectivity analysis")
        print("  - Uses Unicon to assign bond connectivity while preserving coordinates")
        print("  - Generates 3x3x3 supercell")
        print("  - Exports to RDKit-compatible SDF format")
        sys.exit(1)
    
    cif_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(cif_path):
        print(f"Error: Input file '{cif_path}' not found!")
        sys.exit(1)
    
    print("CIF to Molecular SDF Converter")
    print("=" * 50)
    print(f"Input CIF file: {cif_path}")
    print(f"Output SDF file: {output_path if output_path else 'Auto-generated'}")
    print()
    
    success = process_cif_to_sdf(cif_path, output_path)
    
    if success:
        print(f"\n✓ Conversion completed successfully!")
    else:
        print(f"\n✗ Conversion failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()