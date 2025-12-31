#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# scripts/build_pipdb_ref_lengths.py
import argparse, csv, glob, os, sys
try:
    import pysam
except Exception:
    pysam = None

def from_fasta(paths):
    import gzip
    out = {}
    for p in paths:
        op = gzip.open if p.endswith((".gz",".gzip")) else open
        with op(p, "rt", encoding="ascii", errors="ignore") as fh:
            name = None; L = 0
            for line in fh:
                if line.startswith(">"):
                    if name is not None: out[name] = out.get(name, 0) or L
                    name = line[1:].strip().split()[0]; L = 0
                else:
                    L += len(line.strip())
            if name is not None: out[name] = out.get(name, 0) or L
    return out

def from_bam(paths):
    if pysam is None:
        sys.exit("pysam required if building from BAM. pip/conda install pysam.")
    out = {}
    for p in paths:
        with pysam.AlignmentFile(p, "rb") as bam:
            for ref, ln in zip(bam.references, bam.lengths):
                out[ref] = ln
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta_glob", default="")
    ap.add_argument("--bam_glob", default="read_map_to_pipDB_bowtie2/*.sorted.bam")
    ap.add_argument("--out", default="results/tables/pipdb_ref_lengths.tsv")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    refs = {}
    if args.fasta_glob:
        fa = glob.glob(args.fasta_glob)
        if not fa: sys.exit(f"No FASTA matched: {args.fasta_glob}")
        refs = from_fasta(fa)
    else:
        bams = glob.glob(args.bam_glob)
        if not bams: sys.exit(f"No BAMs matched: {args.bam_glob}")
        refs = from_bam(bams)

    with open(args.out, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        for k,v in sorted(refs.items()):
            w.writerow([k, v])
    print(f"Wrote {len(refs)} refs -> {args.out}")

if __name__ == "__main__":
    main()
