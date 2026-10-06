#!/usr/bin/env python3
"""
fetch_cod.py -- Build the COD input tree and its id list.

The retrieval step was originally done by hand before the run, and the id list
(``cod_ids.txt``, ~527k entries) was kept as a data file. This script makes that
step reproducible: it walks a local COD mirror, or drives the COD REST API to
discover entries, and writes both the id list and the sharded directory layout
that ``cif2pdb.py`` expects::

    <root>/<id[0]>/<id[1:3]>/<id[3:5]>/<id>.cif

Two modes
---------
``--from-mirror DIR``
    Index an already-downloaded COD mirror (fast, offline). This is the mode
    used for the reported run.

``--from-api``
    Query the COD REST endpoint for entry ids. Network required.

Usage
-----
    python fetch_cod.py --from-mirror /data/cod --out data/cod_ids.txt
    python fetch_cod.py --from-api --limit 1000 --out data/cod_ids.txt

Note
----
The full COD holds ~527k entries (as of 2025-08). Only the id list is meant to
be committed to the repository; the CIF files themselves are downloaded on
demand.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

COD_API = "https://www.crystallography.net/cod/result"


def iter_ids_from_mirror(root):
    """Yield COD ids found under a local mirror, in sorted order."""
    root = Path(root)
    ids = set()
    for cif in root.rglob("*.cif"):
        stem = cif.stem
        if stem.isdigit():
            ids.add(stem)
    yield from sorted(ids)


def iter_ids_from_api(limit=None, timeout=30):
    """Yield COD ids from the REST API (``format=lst`` returns bare ids)."""
    url = f"{COD_API}?format=lst"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            payload = response.read().decode("utf-8", errors="ignore")
    except (urllib.error.URLError, TimeoutError) as exc:
        print(f"Could not reach the COD API: {exc}", file=sys.stderr)
        return

    count = 0
    for token in payload.replace(",", "\n").split():
        token = token.strip().strip('"[]')
        if not token.isdigit():
            continue
        yield token
        count += 1
        if limit and count >= limit:
            return


def write_id_list(ids, out_path):
    """Write one id per line. Returns the number written."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", encoding="utf-8") as handle:
        for cod_id in ids:
            handle.write(f"{cod_id}\n")
            n += 1
    return n


def plan_downloads(ids, root):
    """Return the (url, destination) pairs for a COD download run."""
    plan = []
    for cod_id in ids:
        dest = (Path(root) / cod_id[0] / cod_id[1:3] / cod_id[3:5]
                / f"{cod_id}.cif")
        url = f"https://www.crystallography.net/cod/{cod_id}.cif"
        plan.append((url, dest))
    return plan


def build_parser():
    ap = argparse.ArgumentParser(
        description="Build the COD id list / input tree for the pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--from-mirror", metavar="DIR",
                        help="index an existing local COD mirror")
    source.add_argument("--from-api", action="store_true",
                        help="query the COD REST API for entry ids")
    ap.add_argument("--out", default="data/cod_ids.txt", help="output id list")
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after N ids (useful for a smoke test)")
    ap.add_argument("--show-plan", action="store_true",
                    help="print the download plan as JSON and exit")
    ap.add_argument("--download-root", default=None,
                    help="with --show-plan, root of the target CIF tree")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.from_mirror:
        ids = list(iter_ids_from_mirror(args.from_mirror))
    else:
        ids = list(iter_ids_from_api(args.limit))

    if args.limit:
        ids = ids[:args.limit]

    if not ids:
        print("No COD ids found.", file=sys.stderr)
        return 1

    if args.show_plan:
        plan = plan_downloads(ids, args.download_root or ".")
        print(json.dumps([{"url": u, "dest": str(d)} for u, d in plan], indent=2))
        return 0

    n = write_id_list(ids, args.out)
    print(f"Wrote {n} COD id(s) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
