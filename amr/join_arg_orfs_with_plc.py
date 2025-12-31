#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

import argparse, csv, os, re, sys

def contig_from_orf(orf_id: str) -> str:
    tok = orf_id.split()[0]                       # keep before first space
    m = re.match(r'(.+)_\d+$', tok)               # strip trailing _<digits> (ORF number)
    return m.group(1) if m else tok

def derive_aa_len(orf_id: str, fallback: str):
    if fallback and fallback.isdigit():
        return int(fallback)
    parts = [p.strip() for p in orf_id.split("#")]
    if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit():
        try:
            return (abs(int(parts[2]) - int(parts[1])) + 1) // 3
        except:
            return ""
    return ""

def main():
    ap = argparse.ArgumentParser(description="Join RGI ORFs with PLC class + contig metrics.")
    ap.add_argument("--rgi", default="results/tables/rgi_orfs.tsv")
    ap.add_argument("--plc", default="results/tables/plc_class_by_contig.tsv")
    ap.add_argument("--metrics", default="results/tables/contig_metrics.tsv")
    ap.add_argument("--out", default="results/tables/arg_orfs_with_plc.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # load PLC
    plc = {}
    with open(args.plc) as fh:
        next(fh)
        for contig, has_rel, has_or, pip, mr, klass in csv.reader(fh, delimiter="\t"):
            plc[contig] = (klass, int(has_rel), int(has_or), int(pip))

    # load metrics
    metrics = {}
    with open(args.metrics) as fh:
        next(fh)
        for contig, L, gc in csv.reader(fh, delimiter="\t"):
            try:
                metrics[contig] = (int(L), float(gc))
            except:
                continue

    kept = 0
    total = 0
    with open(args.out, "w", newline="") as out, open(args.rgi) as rfh:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample","contig_id","orf_id","arg_class","arg_family","orf_length_aa",
                    "plc_class","has_relaxase","has_orit","pipdb_hit","length","gc"])
        next(rfh)
        for sample, contig_raw, orf_id, aro, arg_class, arg_family, model_type, perc_id, orf_len in csv.reader(rfh, delimiter="\t"):
            total += 1
            contig = contig_from_orf(orf_id)
            if contig not in plc or contig not in metrics:
                # fallback: try stripping ORF-like suffix from contig column too
                tok = contig_raw.split()[0]
                m = re.match(r'(.+)_\d+$', tok)
                contig2 = m.group(1) if m else tok
                if contig2 not in plc or contig2 not in metrics:
                    continue
                contig = contig2
            aa = derive_aa_len(orf_id, orf_len)
            L, gc = metrics[contig]
            klass, hr, ho, ph = plc[contig]
            w.writerow([sample, contig, orf_id, arg_class, arg_family, aa, klass, hr, ho, ph, L, gc])
            kept += 1

    print(f"RGI rows in: {total} | joined rows: {kept} | out: {args.out}")

if __name__ == "__main__":
    main()
