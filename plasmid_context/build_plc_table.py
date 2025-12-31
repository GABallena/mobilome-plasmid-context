#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.
#
# Build a per-contig "plasmid-like context" (PLC) table by combining:
# 1) MOB-suite/mobrecon features per contig
# 2) PIPdb contig hit calls (optional)

import argparse, csv, os

def main():
    ap = argparse.ArgumentParser(description="Build PLC class per contig from MOB-suite features and optional PIPdb hits.")
    ap.add_argument("--mobrecon", default="results/tables/plc_features.mobrecon.tsv",
                    help="TSV with header; must include: sample, contig, molecule, mobility, has_relaxase, has_orit, mpf_type")
    ap.add_argument("--pip_hits", default="results/tables/pipdb_contig_hits.tsv",
                    help="Optional TSV from bam_to_pipdb_hits.py (contig_id plus hit metrics).")
    ap.add_argument("--out", default="results/tables/plc_class_by_contig.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    feats = {}  # contig -> dict
    with open(args.mobrecon) as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = {h:i for i,h in enumerate(header)}
        need = ["sample","contig","molecule","mobility","has_relaxase","has_orit","mpf_type"]
        missing = [k for k in need if k not in cols]
        if missing:
            raise SystemExit(f"Missing columns in {args.mobrecon}: {missing}")

        for line in f:
            row = line.rstrip("\n").split("\t")
            contig = row[cols["contig"]]
            molecule = row[cols["molecule"]]
            mobility = row[cols["mobility"]]
            has_rel = row[cols["has_relaxase"]]
            has_or  = row[cols["has_orit"]]

            d = feats.setdefault(contig, {
                "has_relaxase":0, "has_orit":0, "mobrecon_pos":0,
                "pred_mobility":"", "molecule":""
            })
            d["has_relaxase"] |= int(has_rel) if str(has_rel).isdigit() else 0
            d["has_orit"]     |= int(has_or)  if str(has_or).isdigit()  else 0
            d["mobrecon_pos"] |= 1 if str(molecule).strip().lower() == "plasmid" else 0
            if mobility and mobility != "-":
                d["pred_mobility"] = mobility
            d["molecule"] = molecule

    pip = set()
    if args.pip_hits and os.path.exists(args.pip_hits) and os.path.getsize(args.pip_hits) > 0:
        with open(args.pip_hits) as f:
            header = f.readline().rstrip("\n").split("\t")
            cols = {h:i for i,h in enumerate(header)}
            if "contig_id" not in cols:
                raise SystemExit(f"{args.pip_hits} missing contig_id column")
            hit_col = None
            for c in ["pipdb_hit","present","hit"]:
                if c in cols:
                    hit_col = cols[c]; break
            for line in f:
                row = line.rstrip("\n").split("\t")
                contig = row[cols["contig_id"]]
                if hit_col is None:
                    pip.add(contig)
                else:
                    v = row[hit_col].strip()
                    if v in ("1","true","TRUE","True","yes","YES"):
                        pip.add(contig)

    def classify(d, contig):
        has_rel = d.get("has_relaxase",0)
        has_or  = d.get("has_orit",0)
        mobility= (d.get("pred_mobility","") or "").strip().lower()
        piphit  = 1 if contig in pip else 0

        if mobility == "conjugative":
            cls = "conjugative"
        elif mobility == "mobilizable" or (has_rel or has_or):
            cls = "mobilizable"
        elif mobility in ("non-mobilizable","non_mobilizable") or piphit==1 or d.get("mobrecon_pos",0)==1:
            cls = "non_mobilizable"
        else:
            cls = "chromosomal"
        return cls, piphit

    with open(args.out, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["contig_id","has_relaxase","has_orit","pipdb_hit","mobrecon_pos","plc_class"])
        n=0
        for contig,d in sorted(feats.items()):
            cls, piphit = classify(d, contig)
            w.writerow([contig, d.get("has_relaxase",0), d.get("has_orit",0), piphit, d.get("mobrecon_pos",0), cls])
            n+=1
    print(f"Wrote {n} contigs -> {args.out}")

if __name__ == "__main__":
    main()
