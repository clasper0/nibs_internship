#!/usr/bin/env python3
"""
generate_sbatch.py -- Emit the SLURM job scripts for the COD pipeline.

Two kinds of job are generated:

``--mode structure``
    One job per COD entry: CIF -> PDB -> SDF -> split -> Unicon -> assemble.
    Mirrors the ``cif2pdb.sh`` / ``gen_cif2sdf_sbatch.py`` driver.

``--mode dimer``
    One job per (reference fragment, probe fragment) pair: scan the assembled
    SDFs for non-covalent dimers. Mirrors ``gen_cod_sbatch.py``.

The generated jobs are written into the output tree so that the companion
``submit_all.py`` can walk that tree and submit them.

Usage
-----
    # one job per structure
    python generate_sbatch.py --mode structure --cod-ids data/cod_ids.txt

    # one job per fragment pair (314 x 17 = 5338 jobs)
    python generate_sbatch.py --mode dimer --out-root /scratch/2507/output

    # then
    python submit_all.py <out-root>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Default cluster settings, matching the NIBS "honda" partition used for the
# original run. Override with the CLI flags below for a different cluster.
DEFAULT_PARTITION = "honda"
DEFAULT_QOS = "short_many"
DEFAULT_TIME = "200:00"
DEFAULT_CONDA_ENV = "csd"


def structure_job(cod_id, cif_file, sdf_dir, scripts_dir, python_bin="python",
                  partition=DEFAULT_PARTITION, qos=DEFAULT_QOS, time=DEFAULT_TIME,
                  env=DEFAULT_CONDA_ENV):
    """The per-structure COD conversion job."""
    return f"""#!/bin/bash
#SBATCH --job-name={cod_id}
#SBATCH --output={sdf_dir}/{cod_id}.log
#SBATCH --partition={partition}
#SBATCH --ntasks=1
#SBATCH --qos={qos}
#SBATCH --time={time}

set -euo pipefail

source ~/.bashrc
conda activate {env}

date
hostname

# 1. CIF -> symmetry-expanded PDB
{python_bin} {scripts_dir}/cif2pdb.py "{cif_file}" "{sdf_dir}/{cod_id}.pdb"

# 2. PDB -> SDF with bond orders (structconvert, or Jmol as a fallback)
structconvert "{sdf_dir}/{cod_id}.pdb" "{sdf_dir}/{cod_id}.sdf"

# 3. split into one PDB per molecule, perceive bonds with Unicon, reassemble
{python_bin} {scripts_dir}/sdf2pdbs.py "{sdf_dir}/{cod_id}.sdf"
bash {scripts_dir}/unicon-sdf.sh
{python_bin} {scripts_dir}/assemble_cell.py {cod_id} \\
    --cif "{cif_file}" --out "{sdf_dir}/{cod_id}_supercell.sdf"

rm -f "{sdf_dir}/{cod_id}.pdb"
date
"""


def dimer_job(pair_name, pair_dir, scripts_dir, sdf_dir, frag1, frag2,
              python_bin="python", partition=DEFAULT_PARTITION, qos=DEFAULT_QOS,
              time=DEFAULT_TIME, env=DEFAULT_CONDA_ENV):
    """The per-fragment-pair dimer extraction job."""
    return f"""#!/bin/bash
#SBATCH --job-name={pair_name}
#SBATCH --output={pair_dir}/dimer_search.log
#SBATCH --partition={partition}
#SBATCH --ntasks=1
#SBATCH --qos={qos}
#SBATCH --time={time}

set -euo pipefail

source ~/.bashrc
conda activate {env}

date
hostname

{python_bin} {scripts_dir}/dimer_search.py \\
    --frag1 {frag1} --frag2 {frag2} \\
    --sdf-dir "{sdf_dir}" \\
    --out "{pair_dir}"

date
"""


def generate_structure_jobs(cod_ids, cif_root, out_root, scripts_dir,
                            sbatch_dir, **kwargs):
    """Write one sbatch file per COD id. Returns the number written."""
    written = 0
    for cod_id in cod_ids:
        cif_file = cif_root / cod_id[:1] / cod_id[1:3] / cod_id[3:5] / f"{cod_id}.cif"
        sdf_dir = out_root / cod_id[:1] / cod_id[1:3] / cod_id[3:5]
        sdf_dir.mkdir(parents=True, exist_ok=True)

        sbatch_dir.mkdir(parents=True, exist_ok=True)
        job_path = sbatch_dir / f"{cod_id}.sbatch"
        job_path.write_text(structure_job(cod_id, cif_file, sdf_dir, scripts_dir, **kwargs))
        written += 1
    return written


def generate_dimer_jobs(ref_frags, probe_frags, out_root, scripts_dir, sdf_dir, **kwargs):
    """Write one sbatch file per fragment pair. Returns the number written."""
    written = 0
    for ref in ref_frags:
        for probe in probe_frags:
            pair_name = f"{ref}_{probe}"
            pair_dir = out_root / ref / pair_name
            pair_dir.mkdir(parents=True, exist_ok=True)

            job_path = pair_dir / "dimer_search.sbatch"
            job_path.write_text(
                dimer_job(pair_name, pair_dir, scripts_dir, sdf_dir, ref, probe, **kwargs)
            )
            written += 1
    return written


def build_parser():
    ap = argparse.ArgumentParser(
        description="Generate SLURM job scripts for the COD benchmark pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--mode", choices=["structure", "dimer"], required=True)
    ap.add_argument("--cod-ids", default="data/cod_ids.txt",
                    help="file with one COD id per line (structure mode)")
    ap.add_argument("--cif-root", default=str(Path.home() / "2507" / "cif"),
                    help="root of the downloaded CIF tree (structure mode)")
    ap.add_argument("--out-root", default=str(Path.home() / "2507" / "output"),
                    help="root of the output tree")
    ap.add_argument("--sbatch-dir", default="jobs",
                    help="where to write the generated sbatch files (structure mode)")
    ap.add_argument("--sdf-dir", default=str(Path.home() / "2507" / "sdf"),
                    help="directory of assembled SDFs to scan (dimer mode)")
    ap.add_argument("--scripts-dir", default=str(Path(__file__).resolve().parent),
                    help="directory holding the pipeline scripts")
    ap.add_argument("--partition", default=DEFAULT_PARTITION)
    ap.add_argument("--qos", default=DEFAULT_QOS)
    ap.add_argument("--time", default=DEFAULT_TIME)
    ap.add_argument("--env", default=DEFAULT_CONDA_ENV, help="conda environment name")
    ap.add_argument("--python", default="python", help="python interpreter")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    kwargs = dict(partition=args.partition, qos=args.qos, time=args.time,
                  env=args.env, python_bin=args.python)

    if args.mode == "structure":
        cod_ids_path = Path(args.cod_ids)
        if not cod_ids_path.is_file():
            print(f"Error: {cod_ids_path} not found", file=sys.stderr)
            return 1
        cod_ids = [line.strip() for line in cod_ids_path.read_text().splitlines()
                   if line.strip()]
        n = generate_structure_jobs(
            cod_ids, Path(args.cif_root), Path(args.out_root),
            args.scripts_dir, Path(args.sbatch_dir), **kwargs
        )
        print(f"Generated {n} structure job(s) in {args.sbatch_dir}")
        return 0

    # dimer mode needs the fragment library
    sys.path.insert(0, args.scripts_dir)
    from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS

    n = generate_dimer_jobs(
        list(REF_FRAGMENTS), list(PROBE_FRAGMENTS),
        Path(args.out_root), args.scripts_dir, args.sdf_dir, **kwargs
    )
    print(f"Generated {n} dimer job(s) in {args.out_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
