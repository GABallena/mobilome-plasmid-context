#!/usr/bin/env bash
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

set -euo pipefail
cd "$(dirname "$0")/.." 2>/dev/null || cd .

# -------- config (edit if you need) --------
MIN_BREADTH=0.70     # for presence calling if we must derive it
MIN_DEPTH=1.0
N_TOP=9              # how many "top backbones" to list
OUTDIR="${OUTDIR:-results/stats}"
TABLES="${TABLES:-results/tables}"
# ------------------------------------------

mkdir -p "$OUTDIR"

say() { printf "%s\n" "$*"; }
have() { [[ -s "$1" ]]; }

# ---------- PIPdb: presence (reads→PIPdb) ----------
PIP_CALLS="$TABLES/pipdb_presence_calls.tsv"
PIP_COV="$TABLES/pipdb_read_coverage.tsv"

if ! have "$PIP_CALLS"; then
  if have "$PIP_COV"; then
    say "Deriving PIPdb presence from pipdb_read_coverage.tsv …"
    awk -v OFS='\t' -v b="$MIN_BREADTH" -v d="$MIN_DEPTH" '
      BEGIN{print "sample","plasmid_ref","present","mean_depth","covered_fraction","covered_bases","ref_len"}
      NR>1{
        present = ($6+0 >= b && $3+0 >= d) ? 1 : 0
        print $1,$2,present,$3,$6,$4,$5
      }' "$PIP_COV" > "$PIP_CALLS"
  else
    say "WARN: No PIPdb presence or coverage tables; skipping PIPdb block."
    : > "$PIP_CALLS"
  fi
fi

pip_any="NA"
pip_n="0"
pip_top_lines=""
if have "$PIP_CALLS"; then
  pip_n=$(awk -F'\t' 'NR>1{s[$1]=1}END{print length(s)}' "$PIP_CALLS")
  pip_any=$(awk -F'\t' 'NR>1{p[$1] = p[$1] || ($3==1)}END{n=0; for(k in p) if(p[k]) n++; print n}' "$PIP_CALLS")
  # top backbones
  awk -F'\t' 'NR>1 && $3==1{c[$2]++}END{for(k in c) printf "%s\t%d\n",k,c[k]}' "$PIP_CALLS" \
  | sort -k2,2nr -k1,1 \
  | head -n "$N_TOP" \
  | awk -v n="$pip_n" '{printf "%s in %d/%d sites\n",$1,$2,n}' \
  > "$OUTDIR/pipdb_top_backbones.txt"
  pip_top_lines=$(paste -sd' ; ' "$OUTDIR/pipdb_top_backbones.txt")
fi

# ---------- CRISPR spacer links (contigs↔plasmids) ----------
SP="$TABLES/spacer_edges.tsv"
sp_edges="NA"; sp_spacers="NA"; sp_top_hosts=""; sp_top_plasmids=""
if have "$SP"; then
  sp_edges=$(awk -F'\t' 'NR>1{n++}END{print n+0}' "$SP")
  sp_spacers=$(awk -F'\t' 'NR>1{s+=$4}END{print s+0}' "$SP")
  # top 5 hosts (excluding literal "NA")
  awk -F'\t' 'NR>1 && $2!="NA"{h[$2]+=$4}END{for(k in h) printf "%s\t%d\n",k,h[k]}' "$SP" \
  | sort -k2,2nr | head -n 5 \
  | awk '{printf "%s (%s)\n",$1,$2}' > "$OUTDIR/crispr_top_hosts.txt" || true
  sp_top_hosts=$(paste -sd' ; ' "$OUTDIR/crispr_top_hosts.txt" 2>/dev/null || printf "")
  # top 5 plasmids
  awk -F'\t' 'NR>1{p[$3]+=$4}END{for(k in p) printf "%s\t%d\n",k,p[k]}' "$SP" \
  | sort -k2,2nr | head -n 5 \
  | awk '{printf "%s (%s)\n",$1,$2}' > "$OUTDIR/crispr_top_plasmids.txt" || true
  sp_top_plasmids=$(paste -sd' ; ' "$OUTDIR/crispr_top_plasmids.txt" 2>/dev/null || printf "")
fi

# ---------- MOB-suite features on plasmid-like contigs ----------
PLC="$TABLES/plc_class_by_contig.tsv"
relax_line=""
plc_breakdown=""
if have "$PLC"; then
  # breakdown of plc_class
  awk -F'\t' '
    BEGIN{OFS="="}
    NR>1{c[$6]++} END{sep=""
      for(k in c){printf "%s%s%s",sep,k,c[k]; sep="; "}
      print ""
    }' "$PLC" > "$OUTDIR/plc_breakdown.txt" || true
  plc_breakdown=$(cat "$OUTDIR/plc_breakdown.txt" 2>/dev/null || printf "")
  # relaxase / oriT among non-chromosomal
  awk -F'\t' 'NR>1 && $6!="chromosomal"{n++; r+=$2; o+=$3} END{
      if(n>0) printf "%d plasmid-like contigs with MOB-suite calls, relaxases on %d (%.2f%%) and oriT on %d (%.2f%%)\n",
                     n, r, 100*r/n, o, 100*o/n
    }' "$PLC" > "$OUTDIR/relax_orit_summary.txt" || true
  relax_line=$(cat "$OUTDIR/relax_orit_summary.txt" 2>/dev/null || printf "")
fi

# ---------- OPTIONAL: pathogens & contig→PIPdb summaries ----------
PATH_CALLS="$TABLES/pathogen_presence_calls.tsv"
path_line=""; path_top=""
if have "$PATH_CALLS"; then
  path_n=$(awk -F'\t' 'NR>1{s[$1]=1}END{print length(s)}' "$PATH_CALLS")
  path_any=$(awk -F'\t' 'NR>1{p[$1] = p[$1] || ($3==1)}END{n=0; for(k in p) if(p[k]) n++; print n}' "$PATH_CALLS")
  awk -F'\t' 'NR>1 && $3==1{c[$2]++}END{for(k in c) printf "%s\t%d\n",k,c[k]}' "$PATH_CALLS" \
  | sort -k2,2nr -k1,1 | head -n 5 \
  | awk -v n="$path_n" '{printf "%s in %d/%d sites\n",$1,$2,n}' > "$OUTDIR/pathogen_top_refs.txt"
  path_top=$(paste -sd' ; ' "$OUTDIR/pathogen_top_refs.txt")
  path_line=$(printf "Read recruitment to the pathogen panel detected presence in %s/%s samples (breadth≥70%% & depth≥1.0×).\n" "$path_any" "$path_n")
fi

# ---------- Compose copy-pastable lines ----------
AB="$OUTDIR/abstract_recruitment_text.txt"
{
  echo "=== RECRUITMENT ABSTRACT LINES ==="
  if [[ -n "$relax_line" ]]; then
    echo "Plasmid-context breakdown: $plc_breakdown."
    echo "Across $relax_line"
  fi
  if [[ -n "$pip_n" && "$pip_n" != "0" ]]; then
    printf "Read recruitment to PIPdb backbones detected presence in %s/%s samples (breadth≥70%% & depth≥1.0×).\n" "$pip_any" "$pip_n"
    if [[ -s "$OUTDIR/pipdb_top_backbones.txt" ]]; then
      echo -n "Most prevalent backbones: "
      paste -sd' ; ' "$OUTDIR/pipdb_top_backbones.txt"
      echo "."
    fi
  fi
  if [[ -n "$path_line" ]]; then
    echo "$path_line"
    if [[ -n "$path_top" ]]; then
      echo "Most prevalent pathogen refs: $path_top."
    fi
  fi
  if [[ -n "$sp_edges" && "$sp_edges" != "NA" ]]; then
    printf "CRISPR evidence connected hosts and plasmids through %s edges supported by %s unique spacers.\n" "$sp_edges" "$sp_spacers"
    [[ -n "$sp_top_hosts"    ]] && echo "Top spacer-rich hosts: $sp_top_hosts."
    [[ -n "$sp_top_plasmids" ]] && echo "Top spacer-hit plasmids: $sp_top_plasmids."
  fi
  echo "=== /RECRUITMENT ABSTRACT LINES ==="
} > "$AB"

say "Wrote $AB"
ls -lh "$AB" "$OUTDIR"/pipdb_top_backbones.txt "$OUTDIR"/crispr_top_hosts.txt "$OUTDIR"/crispr_top_plasmids.txt 2>/dev/null || true
