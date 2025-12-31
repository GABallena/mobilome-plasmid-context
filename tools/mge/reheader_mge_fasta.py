#!/usr/bin/env python3
import re

inp  = "MGE_db/mge_db.annot.fasta"
outp = "MGE_db/mge_db.short.fasta"
mapt = "MGE_db/mge_db.idmap.tsv"

allowed = re.compile(r'[^A-Za-z0-9_.-]')

i = 1
seen = set()

def flush(hdr, seq, fout, fmap):
    global i
    if not hdr:
        return
    # hdr looks like >GeneID|Category|Class|Function
    header = hdr[1:].strip()
    parts = header.split("|")
    gid   = parts[0]
    cat   = parts[1] if len(parts) > 1 else "MGE"
    acls  = parts[2] if len(parts) > 2 else "Unknown"
    func  = parts[3] if len(parts) > 3 else "NA"

    sid = f"MGE{i:06d}"
    i += 1
    while sid in seen:
        sid = f"MGE{i:06d}"
        i += 1
    seen.add(sid)

    # sanitize just in case
    sid = allowed.sub("_", sid)[:50]
    # Write new header: short id BEFORE first space; put rich metadata AFTER space
    desc = f"orig={gid} cat={cat} class={acls} func={func}"
    fout.write(f">{sid} {desc}\n")
    s = "".join(seq)
    for k in range(0, len(s), 60):
        fout.write(s[k:k+60] + "\n")
    print(f"{sid}\t{gid}\t{cat}\t{acls}\t{func}", file=fmap)


with open(inp) as fin, open(outp, "w") as fout, open(mapt, "w") as fmap:
    print("short_id\torig_id\tcategory\taclass\tfunction", file=fmap)
    hdr = None
    seq = []
    for line in fin:
        if line.startswith(">"):
            flush(hdr, seq, fout, fmap)
            hdr = line.strip()
            seq = []
        else:
            seq.append(line.strip())
    flush(hdr, seq, fout, fmap)

print(f"Wrote {outp} and map {mapt}")
