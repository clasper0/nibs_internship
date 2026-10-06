# The CSD branch, and why it is not in this repository

The project had two halves. The COD half is what this repository implements
end to end. The CSD half is described here for completeness, because the two
share the same post-processing logic and the shared fragment library.

## What the CSD branch did

Per the workflow slides:

```
CSD data
  -> predefined fragment dimers        (17 residues x 314 splinters)
  -> substructure search               (CSD Python API)
  -> preliminary non-covalent dimers
  -> post processing
       - hydrogen correction           (RDKit)
       - aligning fragments            (RDKit)
       - removing redundant molecules  (RMSD)
  -> pool of fragment dimers
```

The search step used the **CSD Python API** (`ccdc` package), which offers
`Molecule.substructure_search` / `CSDEntry` accessors and a `CSDConformer`
generator. The distance criterion and the post-processing steps were identical
to the COD branch:

* keep a pair when any heavy-atom contact is `< 5 Å`
* repair hydrogens against a SMILES-derived template
* superpose all hits onto a reference frame
* drop duplicates below 1e-5 Å RMSD

## Why it is absent

`csd-python-api` is **licence-gated**. It installs from
`https://conda.ccdc.cam.ac.uk`, and every process that imports it requires an
activated licence:

```
$ micromamba run -n csd python -c "import ccdc"
In order to run the software, you must activate it first.
```

The CSD structures themselves are also not redistributable, so neither the
searched entries nor a reproduction recipe that depends on them can be
published.

## How to reinstate it

If you have CSD entitlement:

```bash
conda install -c https://conda.ccdc.cam.ac.uk csd-python-api=3.4.1
# then activate the licence with the CCDC activation tool
```

The post-processing half already exists in this repository and needs no CSD
licence. `src/dimer_search.py` exposes each step as a standalone function:

| Step | Function |
| --- | --- |
| Substructure search | `search()` — RDKit-based equivalent of the CSD API search |
| Fragment extraction | `construct_frag_pair()` |
| Hydrogen correction | `correct_H()` |
| Alignment | `align_sidemol()` |
| De-duplication | `remove_redundant()` |

So a CSD driver would only need to replace `search()` with the CSD API call and
feed the resulting fragment pairs into the same functions — which is precisely
why the two branches produce comparable test sets.

## What was verified locally

The `csd` environment described by the project's `csd.yml` was rebuilt during
this work. `csd-python-api 3.4.1` installs correctly and all of its compiled
dependencies resolve; only the licence check blocks execution. The COD branch
therefore runs on the same environment, minus the licence requirement.
