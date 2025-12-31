#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# scripts/parse_mobrecon_contigs.py
import os, csv, re, glob

def base_contig(s):
    # e.g., PROJECT-02_S2_k141_422940_flag=0_multi=4.0000_len=1225 -> PROJECT-02_S2_k141_422940
    return re.sub(r'_flag=.*$', '', s)

def truthy(x):
    return x and x.strip() != '-' and x.strip().lower() != 'na'

def main():
    in_glob = "mob_results/*/contig_report.txt"
    out_tsv = "results/tables/plc_features.mobrecon.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    kept = 0
    with open(out_tsv, "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["sample","contig_id","molecule_type","predicted_mobility","has_relaxase","has_orit","mpf_type"])
        for path in sorted(glob.glob(in_glob)):
            # sample folder name
            sample = os.path.basename(os.path.dirname(path))
            with open(path, newline="") as f:
                header = f.readline().rstrip("\n").split("\t")
                # We only need early columns; MOB-suite sometimes mangles late header cols.
                # Use fixed indices based on spec (0-based):
                #  1:molecule_type, 4:contig_id, 11:relaxase_type(s), 13:mpf_type, 15:orit_type(s), 17:predicted_mobility
                for line in f:
                    row = line.rstrip("\n").split("\t")
                    if len(row) < 18:
                        continue
                    molecule = row[1]
                    contigid  = base_contig(row[4])
                    relaxase  = row[11]
                    mpf       = row[13]
                    orit      = row[15]
                    mobility  = row[17]
                    has_rel   = 1 if truthy(relaxase) else 0
                    has_or    = 1 if truthy(orit) else 0
                    w.writerow([sample, contigid, molecule, mobility, has_rel, has_or, mpf if truthy(mpf) else ""])
                    kept += 1
    print(f"Wrote {kept} rows to {out_tsv}")

if __name__ == "__main__":
    main()
