#!/usr/bin/env bash
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

set -euo pipefail

MAP=results/tables/pipdb_meta_keymap.tsv
mkdir -p results/tables

awk -v OFS='\t' -F'\t' '
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i
  need=("id plasmid_acc plasmid_name species_name host_rank2 country replicon_type drugclass arg_WHO length_avg")
  # check minimal headers
  if(!( "plasmid_acc" in h) || !("species_name" in h)){
    print "ERROR: required headers missing" > "/dev/stderr"; exit 1
  }
  next
}
{
  id      = $(h["id"])
  acc     = $(h["plasmid_acc"])      # PSCxxxxx
  pname   = ( "plasmid_name"  in h ? $(h["plasmid_name"])  : "-" )
  species = ( "species_name"  in h ? $(h["species_name"])  : "-" )
  host    = ( "host_rank2"    in h ? $(h["host_rank2"])    : "-" )
  country = ( "country"       in h ? $(h["country"])       : "-" )
  rep     = ( "replicon_type" in h ? $(h["replicon_type"]) : "-" )
  drug    = ( "drugclass"     in h ? $(h["drugclass"])     : "-" )
  who     = ( "arg_WHO"       in h ? $(h["arg_WHO"])       : "-" )
  lenavg  = ( "length_avg"    in h ? $(h["length_avg"])    : "-" )

  # emit a single canonical keymap row keyed by the PSC accession
  if(acc != "" && acc != "\\N"){
    print acc, id, pname, species, host, country, rep, drug, who, lenavg
  }
}
' "$META" | awk -F'\t' '!seen[$1]++' > "$MAP"

# quick check
echo -e "key\tPSC_id\tplasmid_name\tspecies\thost\tcountry\treplicon\tdrugclass\tWHO_ARG\tlen_avg" | cat - "$MAP" | head
wc -l "$MAP"
# Rebuild coverage
COV=results/tables/pipdb_read_coverage.tsv
echo -e "sample\tplasmid_ref\tmean_depth\tcovered_bases\tref_len\tcovered_fraction" > "$COV"

ls read_map_to_pipDB_bowtie2/*.sorted.bam \
| parallel -j 6 '
  b={};
  s=$(basename "$b" .sorted.bam);
  samtools coverage -q 0 -Q 0 "$b" \
  | awk -v OFS="\t" -v S="$s" "NR>1{len=\$3-\$2+1; cf=(len>0?\$5/len:0); printf \"%s\t%s\t%.6f\t%d\t%d\t%.6f\n\", S,\$1,\$7,\$5,len,cf}"
' >> "$COV"

# Annotate with PSC meta (key is plasmid_ref = PSC accession)
ANN=results/tables/pipdb_read_coverage_annot.tsv
awk -F'\t' 'FNR==NR{meta[$1]=$0; next}
BEGIN{
  OFS="\t";
  print "sample","plasmid_ref","mean_depth","covered_bases","ref_len","covered_fraction",
        "PSC_id","plasmid_name","species","host_rank2","country","replicon_type",
        "drugclass","WHO_ARG","psc_len_avg"
}
NR>1{
  k=$2;
  if(k in meta){
    split(meta[k],a,"\t");
    print $0, a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10]
  } else {
    print $0, "-","-","-","-","-","-","-","-","-"
  }
}' "$MAP" "$COV" > "$ANN"

# sanity
head "$ANN"
