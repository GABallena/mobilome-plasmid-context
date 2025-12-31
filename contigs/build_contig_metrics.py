#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

import argparse, glob, gzip, os

FASTA_EXTS = (".fa", ".fasta", ".fna")
GZ_EXTS    = (".gz", ".gzip")

def is_fasta_like(path: str) -> bool:
    # Must end with a fasta-like extension (optionally gz)
    base = path
    if path.endswith(GZ_EXTS):
        base = path[:path.rfind(".")]
    return base.endswith(FASTA_EXTS)

def first_byte_is_gt(path: str) -> bool:
    op = gzip.open if path.endswith(GZ_EXTS) else open
    try:
        with op(path, "rb") as fh:
            b = fh.read(1)
        return b == b">"
    except Exception:
        return False

def fa_iter(path: str):
    op = gzip.open if path.endswith(GZ_EXTS) else open
    # Read text as ASCII; ignore weird bytes if any
    with op(path, "rt", encoding="ascii", errors="ignore") as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name:
            yield name, "".join(seq)

def main():
    ap = argparse.ArgumentParser(description="Compute contig length and GC from FASTA files.")
    ap.add_argument("--fasta_globs", nargs="+",
                    default=[
                        "relabeled_contigs/*_contigs.fa",
                        "relabeled_contigs/*_contigs.fasta",
                        "relabeled_contigs/*_contigs.fna",
                        "relabeled_contigs/*_contigs.fa.gz",
                        "relabeled_contigs/*_contigs.fasta.gz",
                        "relabeled_contigs/*_contigs.fna.gz",
                    ],
                    help="One or more glob patterns for FASTA files.")
    ap.add_argument("--out", default="results/tables/contig_metrics.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Expand globs, dedupe, and keep only real fastas whose first byte is '>'
    files = []
    seen = set()
    for pat in args.fasta_globs:
        for p in glob.glob(pat):
            if p in seen: 
                continue
            seen.add(p)
            if is_fasta_like(p) and first_byte_is_gt(p):
                files.append(p)

    n = 0
    with open(args.out, "w") as out:
        out.write("contig_id\tlength\tgc\n")
        for fa in sorted(files):
            for name, seq in fa_iter(fa):
                L = len(seq)
                if L == 0:
                    continue
                gc = sum(1 for c in seq if c in "GgCc") / L
                out.write(f"{name}\t{L}\t{gc:.6f}\n")
                n += 1
    print(f"Wrote {n} contigs to {args.out} from {len(files)} FASTA files.")

if __name__ == "__main__":
    main()
