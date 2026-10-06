# Constructing CSD and COD Test Sets for Transferable Neural Network Potentials

Code for the NIBS summer internship project carried out in **Dr. Niu Huang's
laboratory** (August 2025).

The group's graph-neural-network potential for protein–ligand interaction
energies had been trained by active and transfer learning, but was still
evaluated only on a limited set of PDB-derived dimers in which inorganic
ligands are poorly represented. This project assembles a **more diverse
benchmark set** from two crystallographic databases — the Cambridge Structural
Database (CSD) and the Crystallography Open Database (COD) — covering on the
order of two million valid entries. The pipeline runs end to end: parsing raw
structures, correcting hydrogens with RDKit, aligning and de-duplicating dimer
conformations, and exporting curated structures in a standardized format.

---

## Pipeline

```
                    ┌──────────────────────── CSD branch ─────────────────────────┐
                    │  CSD Python API  ->  predefined fragment dimer search      │
                    │  (requires a CSD licence; see docs/csd_branch.md)          │
                    └────────────────────────────┬───────────────────────────────┘
                                                 │
                                                 ▼
   ┌──────────────────────── COD branch (this repo) ─────────────────────────┐
   │                                                                          │
   │  COD .cif                                                                │
   │     │  src/cif2pdb.py        symmetry expansion -> (super)cell PDB       │
   │     ▼                                                                    │
   │  PDB ── structconvert / Jmol ──► combined SDF   (bond orders)            │
   │     │                                                                    │
   │     │  src/sdf2pdbs.py       split into one molecule per file            │
   │     ▼                                                                    │
   │  per-molecule PDB ── Unicon ──► per-molecule SDF (bonds + charges)       │
   │     │                                                                    │
   │     │  src/assemble_cell.py  merge + rebuild cell + drop bogus bonds     │
   │     ▼                                                                    │
   │  labelled (super)cell SDF                                                │
   │     │                                                                    │
   │     │  src/dimer_search.py   SMARTS search -> H correction ->            │
   │     ▼                        alignment -> RMSD de-duplication            │
   │  <ref>_<probe>.sdf  :  curated non-covalent fragment dimers              │
   └──────────────────────────────────────────────────────────────────────────┘
```

The post-processing half is shared by both branches: `dimer_search.py`
implements the same steps (substructure match → hydrogen correction →
alignment → RMSD de-duplication) that the CSD workflow performs through the CSD
Python API, so the two databases produce directly comparable test sets.

### Fragment library

`src/all_frags.py` holds the 331 predefined fragments:

| Set | Count | Description |
| --- | --- | --- |
| `REF_FRAGMENTS` | **314** | "Splinters dataset" reference fragments (`frag_1` … `frag_314`) |
| `PROBE_FRAGMENTS` | **17** | Predefined residue / ligand fragments (`NMA`, `ACEM`, `HOH`, …) |

A full run therefore covers `314 × 17 = 5338` fragment pairs.

---

## Repository layout

```
.
├── src/                        # importable pipeline library
│   ├── cif_io.py               # self-contained CIF parser (numpy only)
│   ├── cif2pdb.py              # stage 1: CIF -> symmetry-expanded PDB
│   ├── sdf2pdbs.py             # stage 2: split a combined SDF per molecule
│   ├── assemble_cell.py        # stage 3: merge Unicon output, rebuild cell
│   ├── dimer_search.py         # stage 4: build curated non-covalent dimers
│   └── all_frags.py            # 314 reference + 17 probe fragment definitions
├── scripts/                    # cluster orchestration
│   ├── fetch_cod.py            # build the COD id list / input tree
│   ├── generate_sbatch.py      # emit per-structure and per-pair SLURM jobs
│   ├── submit_all.py           # submit a job tree, with resume
│   ├── unicon-sdf.sh           # Unicon bond perception over the scratch dir
│   └── cif2sdf.sh              # Jmol CIF -> SDF fallback (no licence needed)
├── tests/                      # pytest suite (22 tests)
├── data/
│   ├── samples/                # a small CIF used by the tests
│   └── cod_ids.txt             # the 527k-entry COD id list (see below)
├── docs/
│   ├── pipeline.md             # detailed stage-by-stage description
│   └── csd_branch.md           # what the CSD half needed and why it is absent
├── environment.yml             # conda environment (mirrors the `csd` env)
├── requirements.txt            # pip-level dependencies
└── README.md
```

---

## Quick start

### 1. Environment

The original run used a conda environment named `csd` (CSD Python API + RDKit +
pymatgen). A trimmed, reproducible version is provided:

```bash
micromamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate cod-benchmark
```

Only `numpy`, `scipy`, `rdkit` and `pymatgen` are needed for the COD branch. The
CSD branch additionally needs `csd-python-api` **and a valid CSD licence**.

### 2. Run the test suite

```bash
pytest -q
```

### 3. Convert a single structure

A worked example using the bundled sample (an aluminium phosphate framework):

```bash
python src/cif2pdb.py data/samples/1000000.cif out/1000000.pdb
python src/sdf2pdbs.py out/1000000.sdf            # needs RDKit
bash scripts/unicon-sdf.sh                        # needs Schroedinger Unicon
python src/assemble_cell.py 1000000 \
    --cif data/samples/1000000.cif --out out/1000000_supercell.sdf
```

Then extract dimers for one fragment pair:

```bash
python src/dimer_search.py --frag1 frag_1 --frag2 MSM \
    --sdf-dir out --out dimers --cutoff 5.0
```

which writes `dimers/frag_1_MSM.sdf` (curated, de-duplicated) plus
`dimers/frag_1_MSM_uncorrect.sdf` (hits whose hydrogens could not be repaired).

### 4. Run at scale on a SLURM cluster

```bash
# one job per COD entry
python scripts/generate_sbatch.py --mode structure \
    --cod-ids data/cod_ids.txt --cif-root /scratch/cif --out-root /scratch/out
python scripts/submit_all.py jobs --limit 50

# one job per fragment pair (5338 jobs)
python scripts/generate_sbatch.py --mode dimer \
    --out-root /scratch/dimers --sdf-dir /scratch/out
python scripts/submit_all.py /scratch/dimers --pattern dimer_search.sbatch --limit 50
```

Both generators emit plain `sbatch` files, so you can inspect and edit them
before submitting. `submit_all.py` skips any job whose log or result already
exists, which makes an interrupted run resumable.

---

## Method notes

**Symmetry expansion.** `cif_io.py` implements the CIF reader directly rather
than depending on ASE, so stage 1 needs nothing but numpy and runs on a cluster
login node. The 3×3×3 supercell is built by applying the transformation matrix
to the cell, then wrapping the coordinates back into the cell — the wrap is
required, because the downstream molecule splitter relies on that folding to
recognise complete molecules.

**Intact-molecule detection.** A crystal cell is a periodic repeat, so some
molecules are cut by the cell boundary. `get_molecules()` treats the most
frequently occurring connected-component size as the complete molecule and
discards the rest.

**Aromaticity from geometry.** Crystallographic deposits rarely carry reliable
aromatic flags. `assign_aromaticity_by_planarity()` re-derives aromaticity by
fitting each ring to its best plane: a ring whose atoms deviate by less than
0.05 Å is flagged aromatic.

**Hydrogen correction.** `correct_H()` re-caps the fragments cut out of the cell
using an embedded, all-atom template built from the fragment's SMILES, aligning
the heavy-atom core and transferring the template hydrogen positions.

**De-duplication.** Every dimer is superposed onto a common reference frame via
its central fragment, and pairs whose side fragments agree to better than
1 × 10⁻⁵ Å RMSD are collapsed to a single representative.

**Distance criterion.** Two fragments form a dimer when any heavy-atom contact
is shorter than 5 Å (`--cutoff`).

---

## Known limitations

* **Incomplete bond-type and formal-charge assignment** for the molecules of
  interest. Bond perception comes from Unicon / `structconvert` on
  crystallographic coordinates and is not always chemically correct; the
  affected hits are written to `*_uncorrect.sdf` rather than silently dropped.
* **Search efficiency.** The dimer search iterates over all COD entries per
  fragment pair, which is the dominant cost of a full run. Optimising the
  substructure-matching logic (for example by pre-indexing fragments) is the
  natural next step.
* **The CSD branch is not included**, because it requires a licensed CSD Python
  API installation. See `docs/csd_branch.md` for what it did and how to
  reinstate it.

---

## Data

`data/cod_ids.txt` lists the COD entry identifiers used for the run
(527k entries, matching the COD's holdings as of 2025-08). Only the identifier
list is versioned; the CIF files themselves are not, and
`scripts/fetch_cod.py` rebuilds the expected sharded directory tree from either
a local mirror or the COD REST API.

The CSD was queried for ~1.3 million structures (as of 2025-07-20); those
structures are **not redistributable** and are consequently not part of this
repository.

---

## Acknowledgements

Thanks to **Dr. Niu Huang** and all members of the Huang Lab, in particular
**Yulin Zhang** and **Qiancheng Xia** for guidance and for the fragment
definitions, **Wenxuan Sun** and **Shihan Xie** for generous sharing of data
and scripts, **Luyao Ma** for server guidance, and **Xintong Zhang** for
script support.

## License

Released under the MIT License — see [LICENSE](LICENSE).
