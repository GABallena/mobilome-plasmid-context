#!/usr/bin/env python3
import sys, csv

# paths (edit if needed)
fasta_in  = "MGE_db/mge_db.fasta"
annot_tsv = "MGE_db/mge_db_annotations.txt"
fasta_out = "MGE_db/mge_db.annot.fasta"

# load annotations (TSV columns: Gene_ID, Category, Class, Function)
ann = {}
with open(annot_tsv, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    # normalize column names
    cols = {c.lower(): c for c in r.fieldnames}
    for row in r:
        gid = row[cols.get("gene_id","Gene_ID")].strip()
        cat = row[cols.get("category","Category")].strip()
        cl  = row[cols.get("class","Class")].strip()
        fn  = row[cols.get("function","Function")].strip()
        ann[gid] = (cat, cl, fn)

# stream FASTA and rewrite headers
def write_record(hdr, seq, out):
    if not hdr: return
    # original id = first token after '>'
    gid = hdr.split()[0]
    gid = gid[1:] if gid.startswith(">") else gid
    if gid in ann:
        cat, cl, fn = ann[gid]
        new_hdr = f">{gid}|{cat}|{cl}|{fn}"
    else:
        # fall back to original header if no annotation found
        new_hdr = f">{gid}"
    out.write(new_hdr + "\n")
    # wrap sequence
    for i in range(0, len(seq), 60):
        out.write(seq[i:i+60] + "\n")

with open(fasta_in) as fin, open(fasta_out, "w") as fout:
    hdr, seq = "", []
    for line in fin:
        if line.startswith(">"):
            write_record(hdr, "".join(seq), fout)
            hdr, seq = line.strip(), []
        else:
            seq.append(line.strip())
    write_record(hdr, "".join(seq), fout)

print(f"Wrote annotated FASTA -> {fasta_out}")
