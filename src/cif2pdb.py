#!/usr/bin/env python3
"""
cif2pdb.py -- COD CIF -> symmetry-expanded (super)cell PDB.

Stage 1 of the COD branch of the benchmark construction pipeline.

The only third-party requirement is numpy: the CIF reader lives in the local
``cif_io`` module, so this stage needs no ASE / pymatgen on the cluster login
node. This script is what feeds ``structconvert`` / Unicon downstream.

Usage
-----
    python cif2pdb.py input.cif [output.pdb] [-n 3]

Examples
--------
    python cif2pdb.py 1000000.cif
    python cif2pdb.py 1000000.cif 1000000.pdb -n 2
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cif_io import Cell, make_supercell, read_cif, write_pdb  # noqa: E402


def process_cif(cif_path, output_path=None, n=1):
    """Read ``cif_path``, expand it n x n x n and write a PDB.

    Parameters
    ----------
    cif_path : str
        Input CIF file.
    output_path : str, optional
        Output PDB. Defaults to ``<cif stem>_supercell.pdb``.
    n : int
        Supercell multiplier applied along all three axes. ``n=1`` writes the
        symmetry-expanded unit cell only.

    Returns
    -------
    bool
        True on success.
    """
    print(f"Loading CIF file: {cif_path}")
    try:
        atoms = read_cif(cif_path)
    except Exception as exc:                      # noqa: BLE001 - report and fail soft
        print(f"Error loading CIF file: {exc}")
        return False

    print("Successfully loaded structure:")
    print(f"  - Formula: {atoms.get_chemical_formula()}")
    print(f"  - Number of atoms in complete unit cell: {len(atoms)}")
    if atoms.cell is not None and np.shape(atoms.cell) == (3, 3):
        print(f"  - Cell parameters: {Cell(atoms.cell).cellpar()}")

    if output_path is None:
        output_path = Path(cif_path).with_suffix("").name + "_supercell.pdb"

    # ``make_supercell`` is always applied, even for n=1, because it also wraps
    # the symmetry-generated coordinates back into the cell (wrap=True). The
    # downstream molecule splitter relies on that folding.
    print(f"\nGenerating {n}x{n}x{n} supercell...")
    try:
        P = np.array([[n, 0, 0], [0, n, 0], [0, 0, n]], dtype=int)
        supercell = make_supercell(atoms, P)
    except Exception as exc:                      # noqa: BLE001
        print(f"Error creating supercell: {exc}")
        return False
    print("Supercell created successfully:")
    print(f"  - Unit cell atoms: {len(atoms)}")
    print(f"  - Supercell atoms: {len(supercell)}")
    print(f"  - Expected atoms: {len(atoms) * (n ** 3)}")
    print(f"  - Formula: {supercell.get_chemical_formula()}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    print(f"\nExporting to PDB format: {output_path}")
    try:
        write_pdb(output_path, supercell)
    except Exception as exc:                      # noqa: BLE001
        print(f"Error writing PDB file: {exc}")
        return False

    print("Successfully exported to PDB.")
    return True


def build_parser():
    ap = argparse.ArgumentParser(
        description="Convert a COD CIF file into a symmetry-expanded PDB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("input_cif", help="input .cif file")
    ap.add_argument("output_pdb", nargs="?", default=None,
                    help="output .pdb file (default: <stem>_supercell.pdb)")
    ap.add_argument("-n", "--supercell", type=int, default=1,
                    help="expand to n x n x n (default: 1 = unit cell only)")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not os.path.exists(args.input_cif):
        print(f"Error: Input file '{args.input_cif}' not found!")
        return 1

    print("CIF to Supercell PDB Converter (Complete Unit Cell)")
    print("=" * 50)
    print(f"Input CIF file: {args.input_cif}")
    print(f"Output PDB file: {args.output_pdb or 'Auto-generated'}")
    print()

    if process_cif(args.input_cif, args.output_pdb, args.supercell):
        print("\nConversion completed successfully.")
        return 0
    print("\nConversion failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
