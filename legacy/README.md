# Legacy scripts (original, as-run)

These are the scripts exactly as they were during the summer internship, kept
for provenance. **They are not part of the current pipeline** — use `src/` and
`scripts/` instead. They are here so that the derivation of the cleaned-up
versions is auditable.

## Provenance

* `cif2pdb.py`, `cif2mols_sdf.py`, `cif_to_supercell_sdf.py`, `cif2pdb_old.py`
  — earlier attempts at the CIF → structure conversion. All three later files
  embed a near-identical copy of the same CIF parser, which is why
  `src/cif_io.py` now factors that reader out once.
* `cod_workflow.py`, `cod_workflow_1.0.py`, `cod_workflow_0.0.py`, `workflow.py`
  — successive versions of the dimer workflow.
  `cod_workflow_1.0.py` → `src/dimer_search.py`.
  `cod_workflow.py` differs from `cod_workflow_1.0.py` only in that two
  `if i > 500: break` guards are commented out.
* `test.py` — the earliest exploration script (pymatgen `VoronoiNN` + SMARTS
  matching prototypes).
* `getfinalsdf.py` → `src/assemble_cell.py`.
* `pdb2mols.py` — an alternative molecule splitter used before the
  `sdf2pdbs.py` + Unicon route. Superseded but retained.
* `gen_cod_sbatch.py`, `gen_cif2sdf_sbatch.py` → `scripts/generate_sbatch.py`.
  `gen_cod_sbatch.py` is **truncated** in the original (its final `with open`
  block is incomplete), so it could not run as-is; the rewrite fixes that.
* `sbatch_all_cod_workflow.sh`, `sbatch_all_cif2sdf_work.sh` →
  `scripts/submit_all.py`.
* `cif2pdb.sh` — the local serial driver. The original contains a stray
  trailing `python` statement (syntax garbage) which is fixed in the rewrite.
* `cif2sdf.sh`, `cif2mol.jmol`, `test.jmol` — the Jmol-based conversion path,
  retained as `scripts/cif2sdf.sh` for use without a Schroedinger licence.
* `unicon-sdf.sh` → `scripts/unicon-sdf.sh`.
* `scpcif.sh` — a one-off helper that copied CIF shards to the lab server.
  Kept only as a record; it references a retired account.

## Not carried forward

* `test.ipynb` (43 cells, 475 KB) — the development scratchpad for
  `cod_workflow.py`. Every function in it appears at least twice; the final
  versions all live in `src/dimer_search.py`. It was excluded rather than
  cleaned up because it added no behaviour.
* `all_frags.py` was carried forward into `src/` with **one fix**: the
  `frag_263` SMILES was `O=c1cc([O-])co1`, which RDKit cannot kekulize, so
  `construct_template()` silently returned `None` for that fragment and it
  never matched. See the comment in `src/all_frags.py`.
* `micromamba/` (1 GB of conda package cache) and the loose SDF/PDB/`.mol`
  scratch outputs were deleted — they were run artefacts, not source.
