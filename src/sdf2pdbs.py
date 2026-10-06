#!/usr/bin/env python3
"""
sdf2pdbs.py -- Split a multi-molecule SDF into per-molecule PDB files.

Stage 2 of the COD branch. ``structconvert`` (Schroedinger) or Jmol writes one
SDF that still contains every molecule of the (super)cell stacked together.
Unicon needs *one molecule per file* to perceive bond orders and charges, so
this script explodes the combined SDF into ``<cod_id>_1.pdb``, ``<cod_id>_2.pdb``
... inside a scratch directory.

Usage
-----
    python sdf2pdbs.py <path/to/1000000.sdf> [--tmp-dir DIR]

The COD id is taken from the SDF file name (first seven characters), matching
the ``<id[0]>/<id[1:3]>/<id[3:5]>/<id>.sdf`` layout produced by the retrieval
scripts.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

# RDKit is extremely chatty about the partially-bonded crystal fragments we are
# about to hand to Unicon. Those messages are expected here.
RDLogger.DisableLog("rdApp.*")

DEFAULT_TMP_DIR = Path.home() / "2507" / "tmp"


def split_sdf(sdf_file, cod_id, tmp_dir=DEFAULT_TMP_DIR):
    """Split ``sdf_file`` into per-fragment PDBs under ``tmp_dir``.

    Returns the number of fragments written.
    """
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    supplier = Chem.SDMolSupplier(
        str(sdf_file), removeHs=False, sanitize=False, strictParsing=False
    )
    structure = supplier[0]
    if structure is None:
        raise ValueError(f"Could not parse any molecule from {sdf_file}")

    mols = Chem.GetMolFrags(
        structure, asMols=True, sanitizeFrags=Chem.SanitizeFlags.SANITIZE_NONE
    )

    written = 0
    for i, mol in enumerate(mols, start=1):
        Chem.MolToPDBFile(mol, str(tmp_dir / f"{cod_id}_{i}.pdb"))
        written = i

    print(f"Generated {written} pdb files in {tmp_dir}")
    return written


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sdf_file", help="combined SDF produced by structconvert / Jmol")
    ap.add_argument("--tmp-dir", default=str(DEFAULT_TMP_DIR),
                    help=f"scratch directory for the split PDBs (default: {DEFAULT_TMP_DIR})")
    ap.add_argument("--cod-id", default=None,
                    help="override the COD id (default: first 7 chars of the file name)")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    sdf_file = Path(args.sdf_file)
    if not sdf_file.is_file():
        print(f"Error: {sdf_file} not found", file=sys.stderr)
        return 1

    cod_id = args.cod_id or sdf_file.name[:7]
    try:
        split_sdf(sdf_file, cod_id, args.tmp_dir)
    except Exception as exc:                      # noqa: BLE001
        print(f"Error splitting {sdf_file}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
