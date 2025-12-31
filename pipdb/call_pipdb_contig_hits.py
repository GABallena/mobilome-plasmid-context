#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# scripts/call_pipdb_contig_hits.py
import os, csv, glob, pysam

def aln_len_and_nm(alg):
    # approximate aligned length on ref = sum(M,=,X, D) and on query = sum(M,=,X, I)
    ref_len = 0
    for op,ln in alg.cigartuples or []:
        if op in (0,2,7,8):  # M,D,=,X
            ref_len += ln
    nm = alg.get_tag("NM") if alg.has_tag("NM") else None
    return ref_len, nm

def main():
    bam_glob = "read_map_to_pipDB_minimap2/*_plasmidhits.bam"
    out_tsv  = "results/tables/pipdb_contig_hits.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    hits = {}  # contig_id -> 1
    n_bams = 0
    for bam in sorted(glob.glob(bam_glob)):
        if not os.path.exists(bam): continue
        n_bams += 1
        with pysam.AlignmentFile(bam, "rb") as bf:
            for alg in bf.fetch(until_eof=True):
                if alg.is_unmapped: continue
                qname = alg.query_name  # contig id
                ref_len, nm = aln_len_and_nm(alg)
                if ref_len < 500: 
                    continue
                pid_ok = True
                if nm is not None and ref_len > 0:
                    pid = 1 - (nm / ref_len)
                    pid_ok = pid >= 0.90
                if pid_ok:
                    hits[qname] = 1
    with open(out_tsv, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["contig_id","pipdb_hit"])
        for c in sorted(hits):
            w.writerow([c, 1])
    print(f"Scanned {n_bams} BAMs, wrote {len(hits)} contigs to {out_tsv}")

if __name__ == "__main__":
    main()
