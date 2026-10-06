# Pipeline reference

Stage-by-stage description of the COD branch. Every stage is an independently
runnable script, and all of them default to reading/writing plain files so that
they can be chained in a SLURM job or run by hand.

```
COD .cif
   │
   │ 1. src/cif2pdb.py
   ▼
(super)cell PDB
   │
   │ 2. structconvert  (Schroedinger; scripts/cif2sdf.sh gives a Jmol fallback)
   ▼
combined SDF with bond orders
   │
   │ 3. src/sdf2pdbs.py
   ▼
one PDB per molecule in the scratch dir
   │
   │ 4. Unicon  (scripts/unicon-sdf.sh)
   ▼
one SDF per molecule, with bonds + formal charges
   │
   │ 5. src/assemble_cell.py
   ▼
labelled (super)cell SDF
   │
   │ 6. src/dimer_search.py
   ▼
<ref>_<probe>.sdf : curated non-covalent dimers
```

---

## Stage 1 — `cif2pdb.py` (CIF → structure)

Reads the CIF, applies every symmetry operation to generate the **complete unit
cell**, optionally expands to an *n × n × n* supercell, wraps all coordinates
back into the cell and writes a PDB.

* The CIF reader is `src/cif_io.py`, implemented from scratch on numpy. No ASE
  or pymatgen is required, so this stage runs on a cluster login node.
* The wrap step is not optional: because the cell is periodic, a molecule that
  straddles a cell face is only recognisable as a whole after its atoms have
  been folded back into the cell. Stage 5's molecule detection depends on it.

```bash
python src/cif2pdb.py input.cif output.pdb          # unit cell
python src/cif2pdb.py input.cif output.pdb -n 2     # 2x2x2 supercell
```

Verified against the original implementation: the refactored script produces
byte-identical PDB output.

## Stage 2 — bond perception (external)

Crystallographic files carry coordinates but no reliable bond orders.
`structconvert` (Schroedinger) or Jmol writes an SDF that carries connectivity
for the whole cell.

* `scripts/cif2sdf.sh` drives Jmol (`{3 3 3}` expansion, V3000 output) and
  needs no licence. It is slower and slightly less accurate than
  `structconvert`, which is preferred when available.

## Stage 3 — `sdf2pdbs.py` (split)

Unicon works on one molecule at a time, so the combined SDF is exploded using
RDKit's `GetMolFrags` and each fragment is written as
`<cod_id>_<n>.pdb` into the scratch directory.

## Stage 4 — Unicon (external)

`scripts/unicon-sdf.sh` runs `unicon` over every PDB in the scratch directory:

```
unicon -i <in>.pdb -o <out>.sdf -t single -p single -g 1
```

`-g 1` keeps the original crystallographic coordinates instead of re-embedding,
which is essential — the whole point is to preserve the experimentally observed
geometry.

## Stage 5 — `assemble_cell.py` (reassemble)

Merges the per-molecule SDFs back into one structure, rebuilds the cell in the
**true crystallographic frame** (cell parameters are re-read from the source
CIF, not estimated from the coordinate bounding box), and removes bonds longer
than 5 Å, which are artefacts of distance-based perception across periodic
images.

```bash
python src/assemble_cell.py 1000000 \
    --cif data/samples/1000000.cif --out out/1000000_supercell.sdf -n 1
```

On the sample structure this removes 224 spurious long bonds from 1356 atoms.

## Stage 6 — `dimer_search.py` (curate dimers)

The scientific core. For one (reference, probe) fragment pair:

1. **Split into molecules.** `get_molecules()` identifies the connected
   components of the cell. Because the cell is a periodic repeat, the most
   common component size is the complete molecule; everything else is cut by
   the boundary.
2. **Re-derive aromaticity.** `assign_aromaticity_by_planarity()` flags a ring
   as aromatic when its atoms deviate from their best-fit plane by < 0.05 Å.
3. **Substructure search.** `search()` matches both SMARTS patterns, keeps only
   matches that lie inside a complete molecule, requires the two hits to be in
   *different* molecules, and accepts the pair when any heavy-atom contact is
   shorter than the cutoff (default 5 Å). One hit per cell is kept, for speed.
4. **Hydrogen correction.** `correct_H()` re-caps each fragment against an
   embedded all-atom template built from its SMILES.
5. **Alignment.** `align_sidemol()` superposes every dimer onto a reference
   frame using its central fragment (heavy atoms only).
6. **De-duplication.** `remove_redundant()` drops dimers whose side fragment
   coincides with a later one to better than 1e-5 Å RMSD.

```bash
python src/dimer_search.py --frag1 frag_1 --frag2 MSM \
    --sdf-dir out --out dimers --cutoff 5.0
```

Outputs:

* `dimers/frag_1_MSM.sdf` — curated, de-duplicated dimers, each named by its
  source COD id.
* `dimers/frag_1_MSM_uncorrect.sdf` — pairs whose hydrogens could not be
  repaired. These are kept rather than dropped, so the failure rate is visible.

---

## Orchestration

`scripts/generate_sbatch.py` emits the SLURM job files; `scripts/submit_all.py`
walks that tree and submits, skipping anything already finished.

```bash
# structures: one job per COD entry
python scripts/generate_sbatch.py --mode structure \
    --cod-ids data/cod_ids.txt --cif-root /scratch/cif --out-root /scratch/out
python scripts/submit_all.py jobs --limit 50

# dimers: one job per fragment pair  (314 x 17 = 5338)
python scripts/generate_sbatch.py --mode dimer \
    --out-root /scratch/dimers --sdf-dir /scratch/out
python scripts/submit_all.py /scratch/dimers --pattern dimer_search.sbatch --limit 50
```

The original drivers submitted jobs in a fixed order with a hard-coded queue
throttle; `submit_all.py` keeps the throttle (`--limit`) but adds resume
support, which matters because a full run is 5338 jobs.

---

## Data layout

The COD shards its files by identifier, and the pipeline follows that layout:

```
<cif_root>/<id[0]>/<id[1:3]>/<id[3:5]>/<id>.cif
<out_root>/<id[0]>/<id[1:3]>/<id[3:5]>/<id>_supercell.sdf
<out_root>/<ref_frag>/<ref_frag>_<probe_frag>/<ref_frag>_<probe_frag>.sdf
```

`scripts/fetch_cod.py` builds both the id list and this tree, from a local
mirror or from the COD REST API.
