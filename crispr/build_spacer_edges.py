#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# scripts/build_spacer_edges.py
import argparse, csv, os
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="work/results/plasmid_host_assignments.contig.tsv",
                    help="Tab with columns: spacer_id, sample, contig, host_lineage, plasmid_id, ...")
    ap.add_argument("--out", default="results/tables/spacer_edges.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    groups = defaultdict(set)  # (sample, host, plasmid) -> set(spacer_id)
    with open(args.hits) as f:
        header = f.readline().rstrip("\n").split("\t")
        # tolerate column order (look up by name)
        cols = {name:i for i,name in enumerate(header)}
        need = ["spacer_id","sample","host_lineage","plasmid_id"]
        for k in need:
            if k not in cols:
                raise SystemExit(f"Missing column '{k}' in {args.hits}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sid  = parts[cols["spacer_id"]]
            samp = parts[cols["sample"]]
            host = parts[cols["host_lineage"]]
            pid  = parts[cols["plasmid_id"]]
            groups[(samp,host,pid)].add(sid)

    with open(args.out, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["sample","host_lineage","plasmid_id","n_unique_spacers"])
        for (samp,host,pid), spacers in sorted(groups.items()):
            w.writerow([samp, host, pid, len(spacers)])
    print(f"Wrote {len(groups)} edges -> {args.out}")

if __name__ == "__main__":
    main()
