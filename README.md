# mobilome-plasmid-context

Templates for analyzing **mobile genetic elements (MGEs)** and **plasmid context** in metagenomic / genome-resolved datasets.

This repo is geared toward questions like:
- which ARGs co-occur with putative MGEs
- what is the contig / neighborhood context of a marker gene
- how do mobilome signals differ across groups/sites/time

## What you’ll typically find here
- parsers for annotation outputs (MGE markers, plasmid predictors, integrons, etc.)
- contig-neighborhood extraction utilities
- co-occurrence tables and summary metrics
- plotting helpers (context diagrams, heatmaps, summary bars)

## Inputs
Provide locally (do not commit):
- contig annotations (TSV/GFF)
- contig sequences (FASTA)
- mapping/coverage summaries (TSV)
- sample metadata (CSV)


## License
MIT — see `LICENSE`.
