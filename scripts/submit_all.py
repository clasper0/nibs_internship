#!/usr/bin/env python3
"""
submit_all.py -- Walk a tree of generated .sbatch files and submit them.

Replaces the two near-identical bash drivers from the original project
(``sbatch_all_cod_workflow.sh`` / ``sbatch_all_cif2sdf_work.sh``) with a single
script that also supports resuming a partially finished run.

Usage
-----
    python submit_all.py <root> [--pattern dimer_search.sbatch] \\
        [--limit 50] [--user $USER] [--match honda] [--dry-run]

Behaviour
---------
* Discovers every ``*.sbatch`` under ``<root>`` (sorted for reproducibility).
* Skips a job whose output file already exists, so the run can be resumed.
* Optionally throttles submission so that no more than ``--limit`` jobs of the
  same name are queued for the user at any time.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def queued_jobs(user=None, match=None):
    """Return the number of queued/running jobs for ``user`` (0 if unknown)."""
    cmd = ["squeue", "-h"]
    if user:
        cmd += ["-u", user]
    if match:
        cmd += ["-n", match]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return 0
    if out.returncode != 0:
        return 0
    return len([line for line in out.stdout.splitlines() if line.strip()])


def find_jobs(root, pattern):
    """All .sbatch files under ``root`` matching ``pattern``, sorted."""
    return sorted(Path(root).rglob(pattern))


def skip_reason(job_path):
    """Return a reason string when the job looks already done, else None.

    A job is considered done when the log file named by its ``--output``
    directive exists and is non-empty, or when a sibling result SDF exists.
    """
    try:
        text = job_path.read_text(errors="ignore")
    except OSError:
        return None

    for line in text.splitlines():
        if line.startswith("#SBATCH --output="):
            log = Path(line.split("=", 1)[1].strip())
            if log.is_file() and log.stat().st_size > 0:
                return f"log exists ({log.name})"

    for result in job_path.parent.glob("*.sdf"):
        if result.stat().st_size > 0:
            return f"result exists ({result.name})"
    return None


def build_parser():
    ap = argparse.ArgumentParser(
        description="Submit every generated sbatch file under a root directory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("root", help="root directory to walk")
    ap.add_argument("--pattern", default="*.sbatch", help="glob for job files")
    ap.add_argument("--limit", type=int, default=0,
                    help="max concurrent jobs for this user (0 = unlimited)")
    ap.add_argument("--user", default=os.environ.get("USER"),
                    help="cluster user to throttle against")
    ap.add_argument("--match", default=None,
                    help="only count jobs whose name matches (for --limit)")
    ap.add_argument("--no-resume", action="store_true",
                    help="submit everything, ignoring existing logs/results")
    ap.add_argument("--dry-run", action="store_true",
                    help="print what would be submitted without submitting")
    ap.add_argument("--poll", type=int, default=30,
                    help="seconds between queue checks when throttling (default: 30)")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    jobs = find_jobs(args.root, args.pattern)
    if not jobs:
        print(f"No jobs matching {args.pattern} under {args.root}")
        return 1

    print(f"Found {len(jobs)} job(s) under {args.root}")

    submitted = skipped = 0
    for job in jobs:
        if not args.no_resume and (reason := skip_reason(job)):
            print(f"skip {job.name}: {reason}")
            skipped += 1
            continue

        if args.limit:
            while queued_jobs(args.user, args.match) >= args.limit:
                time.sleep(args.poll)

        if args.dry_run:
            print(f"would submit: {job}")
        else:
            result = subprocess.run(["sbatch", str(job)], capture_output=True, text=True)
            if result.returncode == 0:
                print(f"submitted {job.name}: {result.stdout.strip()}")
            else:
                print(f"FAILED {job.name}: {result.stderr.strip()}", file=sys.stderr)
        submitted += 1

    print(f"\nDone. submitted={submitted} skipped={skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
