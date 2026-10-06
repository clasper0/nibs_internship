#!/usr/bin/env python3
"""
dimer_search.py -- Build non-covalent fragment dimers from a (super)cell SDF.

This is the downstream half of the COD branch, and the part that mirrors the
CSD workflow: given a crystal structure that has already been converted to an
SDF with bond orders (see ``assemble_cell.py``), find every pair of predefined
fragments that sits closer than a distance cutoff, repair the hydrogens using a
SMILES-derived template, superpose all hits onto a common reference frame and
drop duplicates.

Pipeline
--------
    SDF (supercell)
      -> substructure search (SMARTS)       Search()
      -> fragment extraction                construct_frag_pair()
      -> hydrogen correction                correct_H()
      -> alignment onto a reference         align_sidemol()
      -> RMSD de-duplication                remove_redundant()
      -> <frag1>_<frag2>.sdf

The two fragments come from ``all_frags``:
  * ``REF_FRAGMENTS``   -- 314 "splinter" fragments taken from the group's
                           existing PDB-derived benchmark
  * ``PROBE_FRAGMENTS`` -- 17 predefined residue/ligand fragments

Usage
-----
    python dimer_search.py --frag1 frag_1 --frag2 MSM \\
        --sdf-dir <dir of per-structure SDFs> --out <dir> [--cutoff 5.0]

Reference
---------
The CSD half of this workflow uses the CSD Python API for the search step;
here the same logic is expressed with RDKit so that the COD half can run
without a CSD licence.
"""

from __future__ import annotations

import argparse
import copy
import os
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from scipy.spatial.distance import cdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Two fragments count as a non-covalent dimer when any heavy-atom contact is
# shorter than this (Angstrom).
DEFAULT_CUTOFF = 5.0

# RMSD below which two aligned dimers are considered the same complex.
REDUNDANCY_RMSD = 1e-5

# Alignment is considered failed above this RMSD.
ALIGNMENT_RMSD_WARN = 0.5


# --------------------------------------------------------------------------- #
# Bond / molecule helpers
# --------------------------------------------------------------------------- #
def simplify_mol(mol):
    """Set every bond to single and clear aromaticity and formal charges.

    Used for the *substructure* copy of a molecule: it makes SMARTS matching
    robust to the imperfect bond orders produced by crystal-structure bond
    perception. Adapted from ``oddt.toolkits.extras.rdkit.fixer``.
    """
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
    return mol


def get_heavy_indices_from_smarts(smarts):
    """Return indices of atoms in ``smarts`` that are explicitly heavy.

    Wildcards (``*``) and explicit hydrogens are skipped: those are the
    attachment points where the fragment is bonded to the rest of the molecule
    and will be re-capped during hydrogen correction.
    """
    query = Chem.MolFromSmarts(smarts)
    if not query:
        print(f"Warning: could not parse SMARTS: {smarts}")
        return []

    heavy = [
        atom.GetIdx()
        for atom in query.GetAtoms()
        if atom.GetSmarts() not in ("*", "[H]", "[#1]")
    ]
    return heavy


def construct_template(smiles):
    """Build an embedded, all-atom template from a SMILES string.

    Atoms are reordered so that each heavy atom is immediately followed by its
    hydrogens; ``correct_H`` relies on that layout.
    """
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    template = Chem.AddHs(template)
    AllChem.EmbedMolecule(template)

    order = []
    for atom in template.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        order.append(atom.GetIdx())
        for neighbour in atom.GetNeighbors():
            if neighbour.GetAtomicNum() == 1:
                order.append(neighbour.GetIdx())
    return Chem.RenumberAtoms(template, order)


def get_bond_length(mol, atom1_idx, atom2_idx):
    """Distance between two atoms of ``mol`` in Angstrom."""
    conf = mol.GetConformer()
    pos1 = np.array(conf.GetAtomPosition(atom1_idx))
    pos2 = np.array(conf.GetAtomPosition(atom2_idx))
    return float(np.linalg.norm(pos2 - pos1))


# --------------------------------------------------------------------------- #
# Cell -> valid molecules
# --------------------------------------------------------------------------- #
def get_molecules(mol):
    """Split a cell into connected components and identify the intact ones.

    A crystal cell is a periodic repeat, so the most frequently occurring
    fragment size corresponds to the complete molecule; anything else is a
    molecule cut by the cell boundary.

    Returns
    -------
    (atom_mol_dict, valid_molecules)
        ``atom_mol_dict[atom_idx] -> mol_id`` and the list of mol ids that are
        considered complete.
    """
    fragments = Chem.GetMolFrags(mol, asMols=False)

    size_counts = {}
    for fragment in fragments:
        size_counts[len(fragment)] = size_counts.get(len(fragment), 0) + 1

    most_common = max(size_counts.values())
    valid_sizes = {size for size, count in size_counts.items() if count == most_common}

    valid_molecules = []
    atom_mol_dict = {}
    for mol_id, fragment in enumerate(fragments):
        if len(fragment) in valid_sizes:
            valid_molecules.append(mol_id)
        for atom_id in fragment:
            atom_mol_dict[atom_id] = mol_id

    return atom_mol_dict, valid_molecules


def assign_aromaticity_by_planarity(mol, atom_mol_dict, valid_molecules, tolerance=0.05):
    """Flag rings as aromatic when they are geometrically planar.

    Crystal structures rarely carry reliable aromatic flags, so aromaticity is
    re-derived from the coordinates: a ring whose atoms deviate from their best
    plane by less than ``tolerance`` Angstrom is treated as aromatic.
    """
    rings = Chem.GetSymmSSSR(mol)
    optimized_molecules = []
    conf = mol.GetConformer()

    for ring in rings:
        if atom_mol_dict[ring[0]] not in valid_molecules:
            continue                                  # ring sits on a cut fragment

        coords = np.array([list(conf.GetAtomPosition(i)) for i in ring])
        centred = coords - coords.mean(axis=0)
        _, _, vh = np.linalg.svd(centred)
        normal = vh[2, :]
        deviation = np.abs(centred @ normal)

        optimized_molecules.append(atom_mol_dict[ring[0]])

        if deviation.max() < tolerance:
            for atom_idx in ring:
                mol.GetAtomWithIdx(atom_idx).SetIsAromatic(True)
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                    bond.SetIsAromatic(True)
                    bond.SetBondType(Chem.BondType.AROMATIC)

    return mol, (optimized_molecules or list(valid_molecules))


def sanitize_lenient(mol):
    """Sanitize with a reduced operation set, tolerating crystal artefacts.

    ``SANITIZE_KEKULIZE`` and ``SANITIZE_PROPERTIES`` are deliberately skipped:
    they are the two steps that most often reject an inorganic or
    organometallic fragment from a crystallographic database. On failure the
    original molecule is returned unchanged.
    """
    operations = (
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS
        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        | Chem.SanitizeFlags.SANITIZE_CLEANUP
        | Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
    )
    backup = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(mol, sanitizeOps=operations)
        return mol
    except ValueError:
        print("Sanitization failed, keeping the unsanitized molecule.")
        return backup


# --------------------------------------------------------------------------- #
# Substructure search
# --------------------------------------------------------------------------- #
def search(cod_id, mol, substruct1, substruct2, atom_mol_dict, valid_molecules,
           cutoff=DEFAULT_CUTOFF):
    """Find one non-covalent fragment pair inside a single complete molecule.

    Only pairs that (a) live in the same complete molecule, (b) sit in two
    *different* molecules of that cell, and (c) have a heavy-atom contact
    shorter than ``cutoff`` are returned. For efficiency only the first hit per
    cell is kept, matching the behaviour of the CSD workflow.
    """
    match1s = mol.GetSubstructMatches(substruct1, useQueryQueryMatches=True)
    match2s = mol.GetSubstructMatches(substruct2, useQueryQueryMatches=True)

    if not (match1s and match2s):
        if match1s and not match2s:
            print(f"{cod_id}: no substructure {Chem.MolToSmiles(substruct2)}")
        elif match2s and not match1s:
            print(f"{cod_id}: no substructure {Chem.MolToSmiles(substruct1)}")
        else:
            print(f"{cod_id}: neither substructure matched")
        return None

    def collect_heavy(matches):
        """Map heavy-atom index sets -> the full match tuple."""
        seen = {}
        for match in matches:
            if not all(atom_mol_dict[match[i]] == atom_mol_dict[match[i + 1]]
                       for i in range(len(match) - 1)):
                continue                              # match spans two molecules
            if atom_mol_dict[match[0]] not in valid_molecules:
                continue                              # cut fragment
            heavy = {idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() != 1}
            if heavy:
                seen[frozenset(heavy)] = match
        return seen

    heavy_1s = collect_heavy(match1s)
    heavy_2s = collect_heavy(match2s)

    conf = mol.GetConformer()
    hits = []
    for indices_1, match1 in heavy_1s.items():
        coords_1 = np.array([list(conf.GetAtomPosition(i)) for i in indices_1])
        for indices_2, match2 in heavy_2s.items():
            if atom_mol_dict[match1[0]] == atom_mol_dict[match2[0]]:
                continue                              # same molecule, not a dimer
            coords_2 = np.array([list(conf.GetAtomPosition(j)) for j in indices_2])
            if cdist(coords_1, coords_2).min() < cutoff:
                hits.append((match1, match2))

    return hits[0] if hits else None


def construct_frag(mol, indices):
    """Extract the sub-molecule of ``mol`` made of ``indices`` (atoms + bonds)."""
    new_mol = Chem.EditableMol(Chem.Mol())
    index_map = {}
    for new_idx, old_idx in enumerate(indices):
        new_mol.AddAtom(mol.GetAtomWithIdx(old_idx))
        index_map[old_idx] = new_idx
    new_mol = new_mol.GetMol()

    old_conf = mol.GetConformer()
    new_conf = Chem.Conformer(len(indices))
    for new_idx, old_idx in enumerate(indices):
        new_conf.SetAtomPosition(new_idx, old_conf.GetAtomPosition(old_idx))
    new_mol.AddConformer(new_conf)

    editable = Chem.EditableMol(new_mol)
    for old_idx in indices:
        for neighbour in mol.GetAtomWithIdx(old_idx).GetNeighbors():
            neighbour_idx = neighbour.GetIdx()
            if neighbour_idx in indices and old_idx < neighbour_idx:
                bond = mol.GetBondBetweenAtoms(old_idx, neighbour_idx)
                editable.AddBond(index_map[old_idx], index_map[neighbour_idx],
                                 bond.GetBondType())
    return editable.GetMol()


def construct_frag_pair(mol, indices1, indices2):
    """Return the two fragments of a dimer as independent molecules."""
    return construct_frag(mol, indices1), construct_frag(mol, indices2)


# --------------------------------------------------------------------------- #
# Hydrogen correction
# --------------------------------------------------------------------------- #
def correct_H(mol, template, heavy_indices):
    """Rebuild the hydrogens of a fragment from a SMILES-derived template.

    Crystallographic databases frequently omit or misplace hydrogens, and the
    fragments cut out of a cell carry dangling valences. Here every non-heavy
    atom of the fragment is turned into a hydrogen and re-attached following
    the template's heavy-atom layout.

    Returns
    -------
    (success, mol)
    """
    candidate = copy.deepcopy(mol)

    # Every atom that is not part of the matched heavy-atom core becomes H
    for i in range(candidate.GetNumAtoms()):
        if i not in heavy_indices:
            candidate.GetAtomWithIdx(i).SetAtomicNum(1)

    conf = candidate.GetConformer()
    num_heavy = len(heavy_indices)
    if template.GetNumAtoms() < num_heavy:
        return False, candidate

    # Superpose the heavy core onto the template core
    heavy_list = sorted(heavy_indices)
    template_heavy = [i for i in range(template.GetNumAtoms())
                      if template.GetAtomWithIdx(i).GetAtomicNum() != 1]
    if len(template_heavy) != num_heavy:
        return False, candidate

    atom_map = list(zip(heavy_list, template_heavy))
    try:
        rmsd, transform = rdMolAlign.GetAlignmentTransform(
            candidate, template, atomMap=atom_map
        )
    except Exception:                              # noqa: BLE001
        return False, candidate

    rotation = transform[:3, :3]
    translation = transform[:3, 3]

    # Place template hydrogens onto the fragment
    template_conf = template.GetConformer()
    for heavy_idx, template_idx in atom_map:
        fragment_atom = candidate.GetAtomWithIdx(heavy_idx)
        template_atom = template.GetAtomWithIdx(template_idx)

        fragment_h = [n.GetIdx() for n in fragment_atom.GetNeighbors()
                      if n.GetAtomicNum() == 1]
        template_h = [n.GetIdx() for n in template_atom.GetNeighbors()
                      if n.GetAtomicNum() == 1]

        for h_frag, h_tmpl in zip(fragment_h, template_h):
            pos = np.array(template_conf.GetAtomPosition(h_tmpl))
            conf.SetAtomPosition(h_frag, (rotation @ pos + translation).tolist())

    return True, candidate


# --------------------------------------------------------------------------- #
# Alignment and de-duplication
# --------------------------------------------------------------------------- #
def align_sidemol(ref_central_mol, qry_mol, mol_idx, num_heavy, rmsd_warn=ALIGNMENT_RMSD_WARN):
    """Superpose the whole dimer using its central fragment as the anchor."""
    qry_mols = Chem.GetMolFrags(qry_mol, asMols=True)
    if not qry_mols:
        return qry_mol

    atom_map = [(i, i) for i in range(num_heavy)]
    try:
        rmsd, transform = rdMolAlign.GetAlignmentTransform(
            qry_mols[0], ref_central_mol, atomMap=atom_map
        )
    except Exception:                              # noqa: BLE001
        print(f"{mol_idx + 1}: alignment failed")
        return qry_mol

    if rmsd > rmsd_warn:
        print(f"{mol_idx + 1} align error: {rmsd}")

    rotation = transform[:3, :3]
    translation = transform[:3, 3]
    conf = qry_mol.GetConformer()
    for atom_idx in range(qry_mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(atom_idx))
        conf.SetAtomPosition(atom_idx, (rotation @ pos + translation).tolist())
    return qry_mol


def remove_redundant(aligned_mols, num_heavy_2, tol=REDUNDANCY_RMSD):
    """Drop dimers whose aligned side fragment coincides with a later one.

    Returns ``(kept_mols, kept_indices)``.
    """
    atom_map = [[(i, i) for i in range(num_heavy_2)]]
    kept = []
    kept_indices = []

    for i, qry_mol in enumerate(aligned_mols):
        redundant = False
        qry_side = Chem.GetMolFrags(qry_mol, asMols=True)[1]
        for j in range(i + 1, len(aligned_mols)):
            ref_side = Chem.GetMolFrags(aligned_mols[j], asMols=True)[1]
            try:
                rmsd = rdMolAlign.CalcRMS(ref_side, qry_side, map=atom_map)
            except Exception:                      # noqa: BLE001
                continue
            if rmsd < tol:
                redundant = True
                break
        if not redundant:
            kept.append(qry_mol)
            kept_indices.append(i)
    return kept, kept_indices


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def process_structure(sdf_path, cod_id, template_1, template_2,
                      substruct1, substruct2, heavy_indices_1, heavy_indices_2,
                      num_heavy_2, cutoff=DEFAULT_CUTOFF):
    """Run the full dimer extraction for one structure.

    Returns ``(frag1_mol, frag2_mol)`` or ``None`` when nothing was found.
    """
    mol = Chem.SDMolSupplier(str(sdf_path), sanitize=False,
                             removeHs=False, strictParsing=False)[0]
    if mol is None:
        print(f"{cod_id}: could not read {sdf_path}")
        return None

    atom_mol_dict, valid_molecules = get_molecules(mol)
    mol, optimized = assign_aromaticity_by_planarity(mol, atom_mol_dict, valid_molecules)
    mol = sanitize_lenient(mol)

    hit = search(cod_id, mol, substruct1, substruct2, atom_mol_dict, optimized, cutoff)
    if hit is None:
        return None

    match1, match2 = hit
    return construct_frag_pair(mol, match1, match2)


def _load_template(fragment, label):
    smarts = fragment.smarts.replace("=", "")     # bond orders are unreliable here
    smiles = fragment.smiles
    template = construct_template(smiles)
    heavy = get_heavy_indices_from_smarts(smarts)
    print(f"{label}: {smarts} -> {smiles} ({len(heavy)} heavy atoms)")
    return template, Chem.MolFromSmarts(smarts), heavy


def build_parser():
    ap = argparse.ArgumentParser(
        description="Extract non-covalent fragment dimers from (super)cell SDFs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--frag1", default="frag_1", help="reference fragment key (REF_FRAGMENTS)")
    ap.add_argument("--frag2", default="MSM", help="probe fragment key (PROBE_FRAGMENTS)")
    ap.add_argument("--sdf-dir", required=True, help="directory of per-structure SDFs")
    ap.add_argument("--out", default=".", help="output directory (default: cwd)")
    ap.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF,
                    help=f"heavy-atom contact cutoff in Angstrom (default: {DEFAULT_CUTOFF})")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.frag1 not in REF_FRAGMENTS:
        print(f"Unknown --frag1 '{args.frag1}'. "
              f"Valid keys look like frag_1 .. frag_314.")
        return 1
    if args.frag2 not in PROBE_FRAGMENTS:
        print(f"Unknown --frag2 '{args.frag2}'. "
              f"Valid keys: {', '.join(sorted(PROBE_FRAGMENTS))}")
        return 1

    pair_name = f"{args.frag1}_{args.frag2}"
    print(f"Starting workflow for fragment pair: {args.frag1} - {args.frag2}")
    print("=" * 50)

    template_1, substruct1, heavy_indices_1 = _load_template(REF_FRAGMENTS[args.frag1], "Fragment 1")
    template_2, substruct2, heavy_indices_2 = _load_template(PROBE_FRAGMENTS[args.frag2], "Fragment 2")
    num_heavy_2 = len([a for a in Chem.MolFromSmiles(PROBE_FRAGMENTS[args.frag2].smiles).GetAtoms()
                       if a.GetSymbol() != "H"])

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    sdf_files = sorted(Path(args.sdf_dir).glob("*.sdf"))
    print(f"\nScanning {len(sdf_files)} structure(s) in {args.sdf_dir}")
    print("=" * 50)

    pairs, names, uncorrected, uncorrected_names = [], [], [], []
    for sdf_path in sdf_files:
        cod_id = sdf_path.stem
        result = process_structure(
            sdf_path, cod_id, template_1, template_2,
            substruct1, substruct2, heavy_indices_1, heavy_indices_2,
            num_heavy_2, args.cutoff,
        )
        if result is None:
            continue
        frag1_mol, frag2_mol = result

        ok1, corrected_1 = correct_H(frag1_mol, template_1, heavy_indices_1)
        ok2, corrected_2 = correct_H(frag2_mol, template_2, heavy_indices_2)
        if not (ok1 and ok2):
            uncorrected.append(Chem.CombineMols(frag1_mol, frag2_mol))
            uncorrected_names.append(cod_id)
            continue
        pairs.append(Chem.CombineMols(corrected_1, corrected_2))
        names.append(cod_id)

    print(f"\n{len(pairs)} dimer(s) passed hydrogen correction, "
          f"{len(uncorrected)} failed")

    if pairs:
        ref_mol = Chem.Mol(pairs[0])
        Chem.SanitizeMol(ref_mol, catchErrors=True)
        ref_frags = Chem.GetMolFrags(ref_mol, asMols=True)
        ref_central = ref_frags[0]

        aligned = [pairs[0]]
        for i, qry in enumerate(pairs[1:], start=1):
            aligned.append(align_sidemol(ref_central, qry, i, num_heavy_2))

        kept, kept_indices = remove_redundant(aligned, num_heavy_2)
        print(f"Removed {len(aligned) - len(kept)} redundant dimer(s), "
              f"{len(kept)} unique remain")
        final_names = [names[i] for i in kept_indices]
    else:
        kept, final_names = [], []

    corrected_file = out_dir / f"{pair_name}.sdf"
    with Chem.SDWriter(str(corrected_file)) as writer:
        for mol, name in zip(kept, final_names):
            mol.SetProp("_Name", name)
            writer.write(mol)
    print(f"Wrote {corrected_file} ({len(kept)} dimers)")

    if uncorrected:
        uncorrected_file = out_dir / f"{pair_name}_uncorrect.sdf"
        with Chem.SDWriter(str(uncorrected_file)) as writer:
            for mol, name in zip(uncorrected, uncorrected_names):
                try:
                    Chem.SanitizeMol(mol)
                except Exception:                  # noqa: BLE001
                    continue
                mol.SetProp("_Name", name)
                writer.write(mol)
        print(f"Wrote {uncorrected_file} ({len(uncorrected)} pairs)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
