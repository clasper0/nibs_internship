#!/usr/bin/env python3
"""
assemble_cell.py -- Merge Unicon output and rebuild the labelled (super)cell SDF.

Stage 3 (final stage) of the COD branch of the benchmark pipeline.

Input : a scratch directory of per-molecule SDF files written by Unicon
        (see ``unicon-sdf.sh``), each carrying perceived bond orders and
        formal charges but *original* crystal coordinates.
Output: one SDF containing the reassembled cell / supercell, plus a report of
        bonds that were dropped as physically unreasonable.

The cell parameters are read back from the source CIF so that the supercell is
built in the *true* crystallographic frame rather than from an estimated box.

Usage
-----
    python assemble_cell.py <cod_id> --cif <path/to/cod_id.cif> \\
        [--sdf-dir DIR] [--out DIR] [-n 1|2|3]

Examples
--------
    python assemble_cell.py 1000000 --cif data/samples/1000000.cif -n 1
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

DEFAULT_SDF_DIR = Path.home() / "2507" / "tmp"

# Bonds longer than this (in Angstrom) cannot be real; they are an artefact of
# naive distance-based perception across periodic images.
MAX_BOND_LENGTH = 5.0


# --------------------------------------------------------------------------- #
# Cell geometry
# --------------------------------------------------------------------------- #
def parse_cif_cell_params(cif_file_path):
    """Return ``[a, b, c, alpha, beta, gamma]`` from a CIF, or ``None``."""
    tags = {
        "_cell_length_a": 0,
        "_cell_length_b": 1,
        "_cell_length_c": 2,
        "_cell_angle_alpha": 3,
        "_cell_angle_beta": 4,
        "_cell_angle_gamma": 5,
    }
    params = [None] * 6
    try:
        with open(cif_file_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                stripped = line.strip()
                for tag, index in tags.items():
                    if stripped.startswith(tag) and params[index] is None:
                        raw = stripped[len(tag):].strip().split("(")[0]
                        try:
                            params[index] = float(raw)
                        except ValueError:
                            pass
    except OSError as exc:
        print(f"Could not read CIF {cif_file_path}: {exc}")
        return None

    if any(value is None for value in params):
        return None
    return params


def cellpar_to_cell(cellpar):
    """Convert ``[a, b, c, alpha, beta, gamma]`` into a 3x3 cell matrix."""
    a, b, c, alpha, beta, gamma = cellpar
    alpha_r, beta_r, gamma_r = np.radians([alpha, beta, gamma])

    cos_alpha, cos_beta, cos_gamma = np.cos([alpha_r, beta_r, gamma_r])
    sin_gamma = np.sin(gamma_r)

    cell = np.zeros((3, 3))
    cell[0, 0] = a
    cell[1, 0] = b * cos_gamma
    cell[1, 1] = b * sin_gamma
    cell[2, 0] = c * cos_beta
    cell[2, 1] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cell[2, 2] = (c * np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2
                              + 2 * cos_alpha * cos_beta * cos_gamma) / sin_gamma)
    return cell


def estimate_cell_params(positions):
    """Fallback: derive a bounding orthorhombic box from the coordinates."""
    mins = positions.min(axis=0)
    maxs = positions.max(axis=0)
    lengths = np.maximum(maxs - mins, 1.0)
    return [float(lengths[0]), float(lengths[1]), float(lengths[2]), 90.0, 90.0, 90.0]


# --------------------------------------------------------------------------- #
# Molecule handling
# --------------------------------------------------------------------------- #
def combine_sdfs(sdf_dir, pattern="*.sdf"):
    """Combine every SDF in ``sdf_dir`` into a single RDKit molecule."""
    sdf_dir = Path(sdf_dir)
    files = sorted(sdf_dir.glob(pattern))

    final_mol = None
    used = 0
    for sdf_file in files:
        try:
            mol = Chem.SDMolSupplier(
                str(sdf_file), sanitize=False, removeHs=False, strictParsing=False
            )[0]
        except OSError:
            continue
        if mol is None:
            continue
        final_mol = mol if final_mol is None else Chem.CombineMols(final_mol, mol)
        used += 1

    if used == 0:
        print(f"No valid sdf files found in {sdf_dir}. Nothing to do.")
        return None
    print(f"Combined {used} SDF file(s) from {sdf_dir}")
    return final_mol


def make_mol_supercell(mol, P=None, cell_params=None, wrap=True):
    """Build a supercell of an RDKit molecule inside its crystallographic cell."""
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("input molecule is empty")

    conf = mol.GetConformer()
    positions = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
         for i in range(mol.GetNumAtoms())]
    )

    if P is None:
        P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    if cell_params is None:
        cell_params = estimate_cell_params(positions)

    cell_matrix = cellpar_to_cell(cell_params)
    supercell_cell = P @ cell_matrix

    # Lattice translations inside the supercell
    diagonals = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    d_points = diagonals @ P
    mins, maxes = d_points.min(axis=0), d_points.max(axis=0) + 1
    grid = np.stack(np.meshgrid(*[np.arange(mins[i], maxes[i]) for i in range(3)],
                               indexing="ij"), axis=-1).reshape(-1, 3)
    frac = grid @ np.linalg.inv(P)
    keep = np.all(frac < 1 - 1e-10, axis=1) & np.all(frac >= -1e-10, axis=1)
    lattice_points = (frac[keep] @ supercell_cell)

    n_copies = len(lattice_points)
    print(f"Supercell {P[0,0]}x{P[1,1]}x{P[2,2]}: "
          f"{mol.GetNumAtoms()} atoms x {n_copies} lattice points "
          f"= {mol.GetNumAtoms() * n_copies}")

    combined = Chem.RWMol()
    for copy_idx, shift in enumerate(lattice_points):
        index_map = {}
        for orig_idx in range(mol.GetNumAtoms()):
            orig_atom = mol.GetAtomWithIdx(orig_idx)
            new_atom = Chem.Atom(orig_atom.GetAtomicNum())
            if orig_atom.GetFormalCharge() != 0:
                new_atom.SetFormalCharge(orig_atom.GetFormalCharge())
            if orig_atom.GetNumRadicalElectrons() != 0:
                new_atom.SetNumRadicalElectrons(orig_atom.GetNumRadicalElectrons())
            index_map[orig_idx] = combined.AddAtom(new_atom)
        for bond in mol.GetBonds():
            combined.AddBond(
                index_map[bond.GetBeginAtomIdx()],
                index_map[bond.GetEndAtomIdx()],
                bond.GetBondType(),
            )
        # Store the translation for each copy so coordinates can be offset later
        if copy_idx == 0:
            all_shifts = [shift]
        else:
            all_shifts.append(shift)

    supercell = combined.GetMol()

    # Coordinates: original positions + lattice translation, optionally wrapped
    new_conf = Chem.Conformer(supercell.GetNumAtoms())
    inv_cell = np.linalg.inv(supercell_cell)
    cursor = 0
    for shift in all_shifts:
        for orig_idx in range(mol.GetNumAtoms()):
            pos = positions[orig_idx] + shift
            if wrap:
                frac_pos = pos @ inv_cell
                frac_pos -= np.floor(frac_pos)
                pos = frac_pos @ supercell_cell
            new_conf.SetAtomPosition(cursor, [float(pos[0]), float(pos[1]), float(pos[2])])
            cursor += 1
    supercell.RemoveAllConformers()
    supercell.AddConformer(new_conf, assignId=True)

    return supercell


def drop_unreasonable_bonds(mol, max_length=MAX_BOND_LENGTH):
    """Remove bonds longer than ``max_length`` (periodic-image artefacts)."""
    conf = mol.GetConformer()
    to_remove = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        pi = np.array(conf.GetAtomPosition(i))
        pj = np.array(conf.GetAtomPosition(j))
        if np.linalg.norm(pi - pj) > max_length:
            to_remove.append((i, j))

    if not to_remove:
        return mol, 0

    rw_mol = Chem.RWMol(mol)
    for i, j in to_remove:
        rw_mol.RemoveBond(i, j)
    return rw_mol.GetMol(), len(to_remove)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    ap = argparse.ArgumentParser(
        description="Assemble the final (super)cell SDF from Unicon output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("cod_id", help="COD identifier, e.g. 1000000")
    ap.add_argument("--cif", required=True, help="source CIF (for the cell parameters)")
    ap.add_argument("--sdf-dir", default=str(DEFAULT_SDF_DIR),
                    help=f"directory of Unicon per-molecule SDFs (default: {DEFAULT_SDF_DIR})")
    ap.add_argument("--out", default=None,
                    help="output SDF path (default: <cod_id>_supercell.sdf in cwd)")
    ap.add_argument("-n", "--supercell", type=int, default=1,
                    help="supercell multiplier along all axes (default: 1)")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    mol = combine_sdfs(args.sdf_dir)
    if mol is None:
        return 1

    cell_params = parse_cif_cell_params(args.cif)
    if cell_params is None:
        print(f"Warning: could not read cell parameters from {args.cif}; "
              f"falling back to an estimated box.")
    else:
        print(f"Cell parameters from CIF: {cell_params}")

    n = args.supercell
    P = np.array([[n, 0, 0], [0, n, 0], [0, 0, n]], dtype=int)
    supercell = make_mol_supercell(mol, P=P, cell_params=cell_params)

    supercell, n_removed = drop_unreasonable_bonds(supercell)
    if n_removed:
        print(f"Removed {n_removed} bond(s) longer than {MAX_BOND_LENGTH} A")

    print(f"Final: {supercell.GetNumAtoms()} atoms, {supercell.GetNumBonds()} bonds")

    out_path = args.out or f"{args.cod_id}_supercell.sdf"
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(out_path)
    supercell.SetProp("_Name", args.cod_id)
    writer.write(supercell)
    writer.close()
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
