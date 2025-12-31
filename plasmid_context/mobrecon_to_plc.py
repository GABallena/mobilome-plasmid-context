#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

import argparse, csv, glob, os, re
from collections import defaultdict

def norm(s: str) -> str:
    return re.sub(r'[^a-z0-9]+','_',s.strip().lower())

def main():
    ap = argparse.ArgumentParser(description="Classify contigs as chromosomal/putative_plasmid/mobilizable.")
    ap.add_argument("--mob_dir", default="mob_results", help="mob_recon outputs root")
    ap.add_argument("--pip_hits", default="results/tables/pipdb_contig_hits.tsv")
    ap.add_argument("--pid_thresh", type=float, default=90.0, help="PIPdb %id threshold")
    ap.add_argument("--cov_thresh", type=float, default=70.0, help="PIPdb reference coverage threshold (%)")
    ap.add_argument("--out", default="results/tables/plc_class_by_contig.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    has_relax, has_orit, mobrecon_pos = defaultdict(bool), defaultdict(bool), defaultdict(bool)

    for f in glob.glob(os.path.join(args.mob_dir, "**/*"), recursive=True):
        if not (f.endswith(".tsv") or f.endswith(".txt")):
            continue
        with open(f, errors="ignore") as fh:
            header = fh.readline()
            if not header or header.count("\t") < 1:
                continue
            hdr = [norm(h) for h in header.rstrip("\n").split("\t")]
            def idx(*cands):
                for c in cands:
                    if c in hdr: return hdr.index(c)
                return None
            ci = idx("contig_id","contig","seq_id","scaffold","name","sequence")
            rlx_i = idx("relaxase_type","relaxase","mob","mob_family")
            ort_i = idx("orit_type","orit")
            pred_i= idx("prediction","type","molecule_type","classification")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                if not row or ci is None or len(row) <= ci: 
                    continue
                cid = row[ci]
                if rlx_i is not None and rlx_i < len(row):
                    v = row[rlx_i].strip().lower()
                    if v and v != "none":
                        has_relax[cid] = True
                if ort_i is not None and ort_i < len(row):
                    v = row[ort_i].strip().lower()
                    if v and v != "none":
                        has_orit[cid] = True
                if pred_i is not None and pred_i < len(row):
                    v = row[pred_i].strip().lower()
                    if "plasmid" in v:
                        mobrecon_pos[cid] = True

    pip_hit = defaultdict(bool)
    if os.path.exists(args.pip_hits):
        with open(args.pip_hits) as fh:
            next(fh)
            for contig, ref, rlen, pident, cov in csv.reader(fh, delimiter="\t"):
                try:
                    pip_hit[contig] = (float(pident) >= args.pid_thresh and float(cov) >= args.cov_thresh)
                except:
                    pip_hit[contig] = False

    allc = sorted(set(has_relax) | set(has_orit) | set(mobrecon_pos) | set(pip_hit))
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["contig_id","has_relaxase","has_orit","pipdb_hit","mobrecon_pos","plc_class"])
        for c in allc:
            hr, ho = int(has_relax[c]), int(has_orit[c])
            ph, mr = int(pip_hit[c]), int(mobrecon_pos[c])
            if hr and ho:
                klass = "mobilizable"
            elif mr or ph:
                klass = "putative_plasmid"
            else:
                klass = "chromosomal"
            w.writerow([c, hr, ho, ph, mr, klass])
    print(f"Wrote {args.out} | n={len(allc)}")

if __name__ == "__main__":
    main()
