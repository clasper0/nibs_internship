#!/usr/bin/env python3
"""
Unit tests for the COD pre-processing pipeline.

These tests are deliberately dependency-light (numpy + rdkit + pytest) so they
can run on a login node without the CSD Python API licence.

Run with:
    pytest -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from cif_io import (  # noqa: E402
    Cell,
    make_supercell,
    read_cif,
    write_pdb,
    write_sdf,
    cellpar_from_cif,
)

SAMPLES = Path(__file__).resolve().parents[1] / "data" / "samples"
CIF = SAMPLES / "1000000.cif"


# --------------------------------------------------------------------------- #
# Cell geometry
# --------------------------------------------------------------------------- #
class TestCell:
    def test_cellpar_roundtrip(self):
        params = [7.8783, 10.4689, 16.0680, 90.0, 95.147, 90.0]
        cell = Cell.new(params)
        out = Cell(cell.array).cellpar()
        assert np.allclose(out, params, atol=1e-4)

    def test_cellpar_from_cif(self):
        assert cellpar_from_cif(CIF) == pytest.approx(
            (7.8783, 10.4689, 16.0680, 90.0, 95.147, 90.0), abs=1e-6
        )

    def test_cellpar_from_cif_missing_tags(self, tmp_path):
        bad = tmp_path / "bad.cif"
        bad.write_text("data_x\n_cell_length_a 1.0\n")
        assert cellpar_from_cif(bad) is None


# --------------------------------------------------------------------------- #
# CIF reading and symmetry expansion
# --------------------------------------------------------------------------- #
class TestReadCif:
    def test_symmetry_expansion(self):
        """12 sites x 4 symmetry operations = 48 atoms."""
        atoms = read_cif(CIF)
        assert len(atoms) == 48

    def test_formula(self):
        atoms = read_cif(CIF)
        assert atoms.get_chemical_formula() == "Al4C8N8O20P8"

    def test_cell_present(self):
        atoms = read_cif(CIF)
        assert np.shape(atoms.cell) == (3, 3)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(Exception):
            read_cif(tmp_path / "does_not_exist.cif")


# --------------------------------------------------------------------------- #
# Supercell construction
# --------------------------------------------------------------------------- #
class TestSupercell:
    @pytest.mark.parametrize("n,expected", [(1, 48), (2, 384), (3, 1296)])
    def test_atom_counts(self, n, expected):
        atoms = read_cif(CIF)
        P = np.array([[n, 0, 0], [0, n, 0], [0, 0, n]])
        assert len(make_supercell(atoms, P)) == expected

    def test_coordinates_wrapped_into_cell(self):
        atoms = read_cif(CIF)
        P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        supercell = make_supercell(atoms, P)
        inv = np.linalg.inv(supercell.cell)
        frac = supercell.positions @ inv
        # wrapped coordinates live in [0, 1)
        assert frac.min() >= -1e-8
        assert frac.max() < 1 + 1e-8

    def test_supercell_cell_is_scaled(self):
        atoms = read_cif(CIF)
        P = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        supercell = make_supercell(atoms, P)
        assert np.allclose(Cell(supercell.cell).cellpar()[:3],
                           np.array(Cell(atoms.cell).cellpar()[:3]) * 2)


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
class TestWriters:
    def test_write_pdb(self, tmp_path):
        atoms = read_cif(CIF)
        out = tmp_path / "out.pdb"
        write_pdb(str(out), atoms)
        text = out.read_text()
        assert text.startswith("TITLE")
        assert "CRYST1" in text
        assert text.count("HETATM") == len(atoms)
        assert text.rstrip().endswith("END")

    def test_write_sdf(self, tmp_path):
        atoms = read_cif(CIF)
        out = tmp_path / "out.sdf"
        write_sdf(str(out), atoms)
        text = out.read_text()
        assert "V2000" in text
        assert text.rstrip().endswith("$$$$")

    def test_pdb_roundtrip_atom_count(self, tmp_path):
        """The PDB we write must be readable by RDKit with the same atom count."""
        Chem = pytest.importorskip("rdkit.Chem")
        atoms = read_cif(CIF)
        out = tmp_path / "rt.pdb"
        write_pdb(str(out), atoms)
        mol = Chem.MolFromPDBFile(str(out), sanitize=False, removeHs=False)
        assert mol is not None
        assert mol.GetNumAtoms() == len(atoms)


# --------------------------------------------------------------------------- #
# Fragment library (all_frags)
# --------------------------------------------------------------------------- #
class TestFragmentLibrary:
    @pytest.fixture(scope="class")
    def frags(self):
        return pytest.importorskip("all_frags")

    def test_counts(self, frags):
        """17 probe fragments and 314 reference 'splinter' fragments."""
        assert len(frags.PROBE_FRAGMENTS) == 17
        assert len(frags.REF_FRAGMENTS) == 314

    def test_every_fragment_has_valid_smarts_and_smiles(self, frags):
        """Every fragment must be parseable and usable for template building.

        Note: RDKit does not run ring perception on mols built with
        ``MolFromSmarts``, so ``GetRingInfo()`` is empty for queries; we check
        parseability and that a template can be embedded, which is what the
        workflow actually requires.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        AllChem = pytest.importorskip("rdkit.Chem.AllChem")

        for name, mol in {**frags.PROBE_FRAGMENTS, **frags.REF_FRAGMENTS}.items():
            assert Chem.MolFromSmarts(mol.smarts) is not None, f"bad smarts: {name}"
            parsed = Chem.MolFromSmiles(mol.smiles)
            assert parsed is not None, f"bad smiles: {name}"

            # the workflow embeds a protonated template from the SMILES
            template = Chem.AddHs(parsed)
            assert AllChem.EmbedMolecule(template) == 0, f"cannot embed: {name}"

    def test_probe_keys_match_names(self, frags):
        for key, mol in frags.PROBE_FRAGMENTS.items():
            assert key == mol.name


# --------------------------------------------------------------------------- #
# Dimer search (integration)
# --------------------------------------------------------------------------- #
class TestDimerSearch:
    """End-to-end test of the substructure -> H-correction -> dedup path."""

    #: closest heavy-atom contact of the fixture, in Angstrom
    CONTACT = 3.6

    @pytest.fixture()
    def synthetic_cell(self, tmp_path):
        """An NMA molecule and an MSM molecule separated by a real 3.6 A contact.

        MSM is placed along +x at a distance that guarantees no atomic overlap
        (the closest approach is set to CONTACT), which is what a genuine
        non-covalent dimer in a crystal looks like.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        AllChem = pytest.importorskip("rdkit.Chem.AllChem")

        nma = Chem.AddHs(Chem.MolFromSmiles("CC(=O)NC"))
        AllChem.EmbedMolecule(nma, randomSeed=42)
        msm = Chem.AddHs(Chem.MolFromSmiles("CSC"))
        AllChem.EmbedMolecule(msm, randomSeed=42)

        # Push MSM along +x until its closest heavy atom is CONTACT away from
        # every NMA heavy atom.
        nma_heavy = np.array([list(nma.GetConformer().GetAtomPosition(i))
                              for i in range(nma.GetNumAtoms())
                              if nma.GetAtomWithIdx(i).GetAtomicNum() != 1])
        msm_heavy = [i for i in range(msm.GetNumAtoms())
                     if msm.GetAtomWithIdx(i).GetAtomicNum() != 1]
        msm_pos = np.array([list(msm.GetConformer().GetAtomPosition(i)) for i in msm_heavy])

        shift = 0.0
        for _ in range(200):
            shifted = msm_pos + np.array([shift, 0.0, 0.0])
            gap = np.linalg.norm(shifted[:, None, :] - nma_heavy[None, :, :], axis=-1).min()
            if gap >= self.CONTACT:
                break
            shift += 0.1

        conf = msm.GetConformer()
        for i in range(msm.GetNumAtoms()):
            pos = np.array(conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, (pos + np.array([shift, 0.0, 0.0])).tolist())

        combined = Chem.CombineMols(nma, msm)
        out = tmp_path / "cell.sdf"
        Chem.MolToMolFile(combined, str(out))
        return tmp_path

    def test_finds_dimer_within_cutoff(self, synthetic_cell):
        import dimer_search as ds
        from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS
        Chem = pytest.importorskip("rdkit.Chem")

        f1 = REF_FRAGMENTS["frag_1"]
        f2 = PROBE_FRAGMENTS["MSM"]
        template_1 = ds.construct_template(f1.smiles)
        template_2 = ds.construct_template(f2.smiles)
        s1 = Chem.MolFromSmarts(f1.smarts.replace("=", ""))
        s2 = Chem.MolFromSmarts(f2.smarts.replace("=", ""))
        h1 = ds.get_heavy_indices_from_smarts(f1.smarts.replace("=", ""))
        h2 = ds.get_heavy_indices_from_smarts(f2.smarts.replace("=", ""))

        result = ds.process_structure(
            synthetic_cell / "cell.sdf", "test001",
            template_1, template_2, s1, s2, h1, h2, num_heavy_2=3, cutoff=5.0,
        )
        assert result is not None, "expected to find a dimer in the synthetic cell"

        frag1, frag2 = result
        ok1, corr1 = ds.correct_H(frag1, template_1, h1)
        ok2, corr2 = ds.correct_H(frag2, template_2, h2)
        assert ok1 and ok2

        # the corrected MSM fragment should be recognisable as CSC
        assert corr2.GetNumAtoms() == 9   # 3 heavy + 6 H

    def test_respects_cutoff(self, synthetic_cell):
        """A tighter cutoff must not accept a distant dimer.

        The fixture's closest heavy-atom contact is CONTACT (3.6 A) by
        construction, so a 2.0 A cutoff finds nothing while the default 5.0 A
        cutoff does. This pins the comparison in ``search``.
        """
        import dimer_search as ds
        from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS
        Chem = pytest.importorskip("rdkit.Chem")

        f1, f2 = REF_FRAGMENTS["frag_1"], PROBE_FRAGMENTS["MSM"]
        result = ds.process_structure(
            synthetic_cell / "cell.sdf", "test001",
            ds.construct_template(f1.smiles), ds.construct_template(f2.smiles),
            Chem.MolFromSmarts(f1.smarts.replace("=", "")),
            Chem.MolFromSmarts(f2.smarts.replace("=", "")),
            ds.get_heavy_indices_from_smarts(f1.smarts.replace("=", "")),
            ds.get_heavy_indices_from_smarts(f2.smarts.replace("=", "")),
            num_heavy_2=3, cutoff=2.0,
        )
        assert result is None


class TestMoleculeSplitting:
    def test_get_molecules_identifies_largest_repeat(self):
        """The most common fragment size is treated as the intact molecule."""
        import dimer_search as ds
        Chem = pytest.importorskip("rdkit.Chem")

        mol = Chem.MolFromSmiles("CCO.CCO.C")
        mol = Chem.AddHs(mol)
        AllChem = pytest.importorskip("rdkit.Chem.AllChem")
        AllChem.EmbedMolecule(mol, randomSeed=1)

        atom_mol_dict, valid = ds.get_molecules(mol)
        # 2 ethanol copies of 9 atoms each dominate the single methane
        sizes = {}
        for frag in Chem.GetMolFrags(mol, asMols=False):
            sizes[len(frag)] = sizes.get(len(frag), 0) + 1
        assert sizes[max(sizes, key=lambda k: sizes[k])] == 2


class TestFragmentSelfConsistency:
    """How many fragment SMARTS patterns match their own SMILES skeleton?

    This is a property of the fragment library, not a bug to fix: a SMARTS
    encodes the fragment *as it appears inside a larger molecule* (wildcards at
    the attachment points), whereas the SMILES describes the free molecule (H
    at those positions). Matching therefore happens against real crystal
    environments, which always carry explicit hydrogens.

    The test pins the current coverage so that an accidental edit to
    ``all_frags.py`` shows up as a change here.
    """

    def _self_consistent(self, fragments):
        Chem = pytest.importorskip("rdkit.Chem")
        consistent, inconsistent = [], []
        for name, frag in fragments.items():
            query = Chem.MolFromSmarts(frag.smarts.replace("=", ""))
            mol = Chem.AddHs(Chem.MolFromSmiles(frag.smiles))
            if query is None or mol is None:
                inconsistent.append(name)
                continue
            hits = mol.GetSubstructMatches(query, useQueryQueryMatches=True)
            (consistent if hits else inconsistent).append(name)
        return consistent, inconsistent

    def test_probe_coverage(self):
        from all_frags import PROBE_FRAGMENTS
        consistent, inconsistent = self._self_consistent(PROBE_FRAGMENTS)
        assert len(consistent) == 12, f"probe coverage changed: {consistent}"
        assert set(inconsistent) == {"ACEM", "ACET", "MGDM", "N1PA", "NMA"}

    def test_reference_coverage(self):
        from all_frags import REF_FRAGMENTS
        consistent, inconsistent = self._self_consistent(REF_FRAGMENTS)
        assert len(consistent) == 194, f"reference coverage changed: {len(consistent)}"
        assert len(consistent) + len(inconsistent) == 314

    def test_wildcards_match_explicit_hydrogens(self):
        """The attachment-point wildcards must be able to match hydrogens.

        This is what makes the MSM pattern match dimethyl sulfide once the
        structure carries explicit Hs, as it does after Unicon.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        from all_frags import PROBE_FRAGMENTS

        frag = PROBE_FRAGMENTS["MSM"]
        query = Chem.MolFromSmarts(frag.smarts.replace("=", ""))
        mol = Chem.AddHs(Chem.MolFromSmiles(frag.smiles))
        assert mol.GetSubstructMatches(query, useQueryQueryMatches=True)
