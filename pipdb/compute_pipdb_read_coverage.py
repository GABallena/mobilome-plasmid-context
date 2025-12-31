#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# scripts/compute_pipdb_read_coverage.py
import argparse, csv, glob, os
import pysam

def pileup_ref(bam, ref, length, min_mapq=30):
    covered = 0
    depth_sum = 0
    for col in bam.pileup(ref, 0, length, truncate=True, stepper="samtools",
                          min_mapping_quality=min_mapq, ignore_overlaps=False,
                          ignore_orphans=False, max_depth=1000000):
        d = col.nsegments
        depth_sum += d
        if d > 0: covered += 1
    mean_depth = depth_sum / max(length,1)
    covered_fraction = covered / max(length,1)
    return mean_depth, covered, covered_fraction

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam_glob", default="read_map_to_pipDB_bowtie2/*.sorted.bam")
    ap.add_argument("--ref_lengths", default="results/tables/pipdb_ref_lengths.tsv")
    ap.add_argument("--min_mapq", type=int, default=30)
    ap.add_argument("--out", default="results/tables/pipdb_read_coverage.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ref lengths
    ref_len = {}
    with open(args.ref_lengths) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts)!=2: continue
            ref_len[parts[0]] = int(parts[1])

    with open(args.out, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["sample","plasmid_ref","mean_depth","covered_bases","ref_len","covered_fraction"])
        for bam_path in sorted(glob.glob(args.bam_glob)):
            sample = os.path.basename(bam_path).replace(".sorted.bam","")
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                # iterate only refs present in this BAM header
                for ref, length in zip(bam.references, bam.lengths):
                    L = ref_len.get(ref, length)
                    md, covb, frac = pileup_ref(bam, ref, L, args.min_mapq)
                    w.writerow([sample, ref, f"{md:.6f}", covb, L, f"{frac:.6f}"])
    print(f"Done -> {args.out}")

if __name__ == "__main__":
    main()
