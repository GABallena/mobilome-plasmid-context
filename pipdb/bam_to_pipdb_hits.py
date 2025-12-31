#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

import argparse, collections, csv, glob, os, sys
try:
    import pysam
except Exception as e:
    print("ERROR: pysam is required (conda/pip install pysam).", file=sys.stderr)
    raise

Stats = collections.namedtuple("Stats", "nm_sum denom_sum intervals")

def add_interval(L, s, e):
    # merge [s,e) into sorted, non-overlapping list L (in place)
    if s >= e: 
        return
    out = []
    placed = False
    for a,b in L:
        if e < a and not placed:
            out.append([s,e]); placed = True
        if b < s or e < a:
            out.append([a,b])
        else:
            s = min(s,a); e = max(e,b)
    if not placed:
        out.append([s,e])
    out.sort()
    L[:] = out

def main():
    ap = argparse.ArgumentParser(description="Summarize contig→PIPdb BAMs into best hits.")
    ap.add_argument("--bam_glob", default="read_map_to_pipDB_minimap2/*_plasmidhits.bam")
    ap.add_argument("--mapq", type=int, default=1, help="min MAPQ")
    ap.add_argument("--out", default="results/tables/pipdb_contig_hits.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # load reference lengths from first BAM
    ref_len = {}
    bams = sorted(glob.glob(args.bam_glob))
    if not bams:
        print("No BAM files matched.", file=sys.stderr); sys.exit(1)
    with pysam.AlignmentFile(bams[0], "rb") as f0:
        ref_len = dict(zip(f0.references, f0.lengths))

    per = collections.defaultdict(lambda: Stats(0,0,[]))  # keyed by (contig, ref)
    n_reads = 0

    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as f:
            for r in f.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary: 
                    continue
                if r.mapping_quality < args.mapq:
                    continue
                q = r.query_name
                t = f.get_reference_name(r.reference_id)
                nm = r.get_tag("NM") if r.has_tag("NM") else 0

                denom = 0
                if r.cigartuples:
                    for op, ln in r.cigartuples:
                        if op in (0,7,8,1):    # M, =, X, I (consume query)
                            denom += ln
                        elif op == 2:         # D (consume ref)
                            denom += ln
                s = r.reference_start
                e = r.reference_end

                st = per[(q,t)]
                if denom > 0:
                    if not st.intervals:
                        per[(q,t)] = Stats(nm, denom, [])
                        st = per[(q,t)]
                    else:
                        per[(q,t)] = Stats(st.nm_sum + nm, st.denom_sum + denom, st.intervals)
                        st = per[(q,t)]
                    add_interval(st.intervals, s, e)
                    n_reads += 1

    # choose best target per contig (highest %id; tie → higher coverage)
    best = {}
    for (q,t), st in per.items():
        pid = 100.0 * (st.denom_sum - st.nm_sum) / st.denom_sum if st.denom_sum > 0 else 0.0
        L = ref_len.get(t, 0)
        cov = 0.0
        if L > 0 and st.intervals:
            cov = 100.0 * sum(e - s for s, e in st.intervals) / L
        cur = best.get(q)
        if cur is None or pid > cur[2] or (abs(pid - cur[2]) < 1e-9 and cov > cur[3]):
            best[q] = (t, L, pid, cov)

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["contig_id","ref_id","ref_len","pident","cov_pct"])
        for q, (t, L, pid, cov) in best.items():
            w.writerow([q, t, L, f"{pid:.3f}", f"{cov:.3f}"])

    print(f"BAMs: {len(bams)} | mapped records used: {n_reads} | contigs with hits: {len(best)}")
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
