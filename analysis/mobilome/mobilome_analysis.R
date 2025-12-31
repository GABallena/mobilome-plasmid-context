#!/usr/bin/env Rscript
# Portfolio-safe script (paths/identifiers generalized; inputs not included).

# When run via `Rscript`, set working directory to the script's directory.
try({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if(length(file_arg) == 1){
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}, silent = TRUE)
cat("Working directory:", getwd(), "\n")

# Mobilome analysis script integrating sample metadata, MOB-suite outputs, plasmid clustering (CD-HIT 90%), and Mash MinHash clustering.
# Author: Auto-generated
# Date: 2025-09-06

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(forcats)
  library(ggrepel)

# ---- Config ----
MOBTYPER_TSV <- Sys.getenv("MOBTYPER_TSV", "results/tables/mobtyper_summary.tsv")
PLASMID_MAP_TSV <- Sys.getenv("PLASMID_MAP_TSV", "results/tables/plasmid_contig_map.tsv")
SAMPLE_META_TSV <- Sys.getenv("SAMPLE_META_TSV", "metadata/sample_metadata.tsv")
OUT_DIR <- Sys.getenv("OUT_DIR", "results/figures/mobilome")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

})

# Global bindings for NSE to avoid NOTES / visible binding warnings
utils::globalVariables(c(
  'blank','sample_code','`SAMPLE CODE`','Sample','relaxase_type(s)','rep_type(s)',
  'predicted_mobility','has_relaxase','has_replicon','n_conjugative','n_mobilizable',
  'n_plasmid_like','sample','cdhit_cluster','mash_cluster','length','size_bin'
))

# -----------------------------
# Helper functions
# -----------------------------
read_metadata <- function(path) {
  # Metadata has leading blank column; sanitize names
  md <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  names(md)[1] <- "blank" # discard later
  md <- md %>% select(-blank)
  # Normalize sample code column name variations
  if ("SAMPLE CODE" %in% names(md)) {
    md <- md %>% rename(sample_code = `SAMPLE CODE`)
  } else if ("Sample" %in% names(md)) {
    md <- md %>% rename(sample_code = Sample)
  }
  md <- md %>% mutate(sample_code = str_trim(sample_code))
  md
}

read_all_mobtyper <- function(path) {
  df <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  # Ensure key columns exist
  stopifnot(all(c("sample","sample_id","predicted_mobility","rep_type(s)","relaxase_type(s)") %in% names(df)))
  df
}

# Parse sequence|sample pattern for seq_to_cluster mapping
read_plasmid_seq_clusters <- function(path) {
  df <- read_tsv(path, col_names = c("seq","cdhit_cluster"), show_col_types = FALSE)
  df <- df %>% separate(seq, into = c("sample","plasmid_seq"), sep = "\\|", remove = FALSE)
  df
}

read_mash_assignments <- function(path) {
  # Two columns: cluster_id (numeric) and contig descriptor
  read_tsv(path, col_names = c("mash_cluster","contig"), show_col_types = FALSE)
}

# Extract sample from contig like PROJECT-02_S2_k141_100439_flag=0_multi=20.1837_len=5765
extract_sample_from_contig <- function(contig) {
  str_extract(contig, "^PROJECT-\\d{2}_S\\d{1,2}")
}

# Length from pattern len=#### at end
extract_len <- function(contig) {
  as.numeric(str_match(contig, "len=([0-9]+)$")[,2])
}

# Coverage / multi value from multi= part
extract_multi <- function(contig) {
  as.numeric(str_match(contig, "multi=([0-9.]+)")[,2])
}

summarize_mob_elements <- function(mob_df) {
  mob_df %>%
    mutate(has_relaxase = `relaxase_type(s)` != '-' & !is.na(`relaxase_type(s)`),
           has_replicon = `rep_type(s)` != '-' & !is.na(`rep_type(s)`)) %>%
    group_by(sample) %>%
    summarise(
      n_plasmid_like = n(),
      n_conjugative = sum(predicted_mobility == 'conjugative', na.rm=TRUE),
      n_mobilizable = sum(predicted_mobility == 'mobilizable', na.rm=TRUE),
      n_non_mobilizable = sum(predicted_mobility == 'non-mobilizable', na.rm=TRUE),
      n_relaxase = sum(has_relaxase),
      n_replicon = sum(has_replicon)
    ) %>%
    ungroup() %>%
    mutate(frac_conjugative = n_conjugative / n_plasmid_like,
           frac_mobilizable = n_mobilizable / n_plasmid_like)
}

# -----------------------------
# Paths (adjust if needed)
# -----------------------------
metadata_path <- "PROJECT DATA MASTERLIST - PROJECT Y1 DATA.tsv"
all_mobtyper_path <- "mob_agg/all_mobtyper_with_len.tsv"
seq_to_cluster_path <- "plasmid_clusters/seq_to_cluster.tsv"
mash_assign_path <- "mash_work/cluster_assignments_d010.tsv"

# -----------------------------
# Load data
# -----------------------------
message("Reading metadata ...")
md <- read_metadata(metadata_path)

message("Reading MOB-typer aggregated results ...")
mob_all <- read_all_mobtyper(all_mobtyper_path)

message("Reading plasmid CD-HIT clustering ...")
pl_cluster <- read_plasmid_seq_clusters(seq_to_cluster_path)

message("Reading Mash MinHash clustering ...")
mash_assign <- read_mash_assignments(mash_assign_path) %>%
  mutate(sample = extract_sample_from_contig(contig),
         contig_len = extract_len(contig),
         multi = extract_multi(contig))

# -----------------------------
# Summaries
# -----------------------------
message("Summarizing MOB-suite elements per sample ...")
mob_summary <- summarize_mob_elements(mob_all)

# Plasmid CD-HIT cluster richness per sample
pl_cluster_summary <- pl_cluster %>%
  distinct(sample, cdhit_cluster) %>%
  count(sample, name = "n_cdhit_clusters")

# Mash cluster richness per sample
mash_cluster_summary <- mash_assign %>%
  distinct(sample, mash_cluster) %>%
  count(sample, name = "n_mash_clusters")

# Combine all per-sample metrics
metrics <- mob_summary %>%
  full_join(pl_cluster_summary, by = "sample") %>%
  full_join(mash_cluster_summary, by = "sample") %>%
  arrange(sample)

# Integrate metadata (match sample code SAMPLE-XX to sample like SAMPLE-XX_SYY via prefix)
metrics <- metrics %>% mutate(sample_code = str_extract(sample, "PROJECT-\\d{2}"))

if ("sample_code" %in% names(md)) {
  # Find columns that look like longitude / latitude even if they have padding spaces
  lon_col <- grep('LONGITUDE', names(md), value = TRUE, ignore.case = TRUE)[1]
  lat_col <- grep('LATITUDE', names(md), value = TRUE, ignore.case = TRUE)[1]
  desc_col <- grep('SAMPLE DESCRIPTION', names(md), value = TRUE, ignore.case = TRUE)[1]
  type_col <- grep('SAMPLE TYPE', names(md), value = TRUE, ignore.case = TRUE)[1]
  md_slim <- md %>% select(sample_code,
                           `SAMPLE DESCRIPTION` = all_of(desc_col),
                           `SAMPLE TYPE` = all_of(type_col),
                           LATITUDE = all_of(lat_col),
                           LONGITUDE = all_of(lon_col))
  metrics <- metrics %>% left_join(md_slim, by = "sample_code")
}

# -----------------------------
# Output tables
# -----------------------------
outfile_metrics <- "mobilome_sample_metrics.tsv"
write_tsv(metrics, outfile_metrics)
message(sprintf("Wrote per-sample metrics to %s", outfile_metrics))

# Long-format mobility composition for stacked bar plot
mob_comp <- mob_all %>%
  select(sample, predicted_mobility) %>%
  filter(predicted_mobility %in% c('conjugative','mobilizable','non-mobilizable')) %>%
  count(sample, predicted_mobility, name = "n") %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# Plot mobility composition
p_mobility <- ggplot(mob_comp, aes(x = sample, y = frac, fill = predicted_mobility)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c(conjugative = '#1b9e77', mobilizable = '#d95f02', `non-mobilizable` = '#7570b3')) +
  labs(title = 'Plasmid mobility composition', x = 'Sample', y = 'Fraction', fill = 'Mobility') +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("mobilome_mobility_composition.png", p_mobility, width = 10, height = 4, dpi = 150)

# Relationship between CD-HIT and Mash cluster richness
p_cluster_corr <- metrics %>%
  ggplot(aes(x = n_cdhit_clusters, y = n_mash_clusters, label = sample)) +
  geom_point(color = '#2c7fb8') +
  ggrepel::geom_text_repel(size = 2, max.overlaps = 30) +
  theme_minimal() +
  labs(title = 'Cluster richness: CD-HIT vs Mash', x = 'Distinct CD-HIT clusters', y = 'Distinct Mash clusters')

ggsave("mobilome_cluster_richness_scatter.png", p_cluster_corr, width = 5, height = 4, dpi = 150)

# Mobility vs cluster richness
p_mobility_vs_cdhit <- metrics %>%
  ggplot(aes(x = n_cdhit_clusters, y = frac_conjugative)) +
  geom_point(color = '#1b9e77') +
  theme_bw() +
  labs(title = 'Conjugative fraction vs CD-HIT cluster richness', x = 'CD-HIT clusters', y = 'Fraction conjugative')

ggsave("mobilome_conjugative_vs_cdhit.png", p_mobility_vs_cdhit, width = 5, height = 4, dpi = 150)

# -----------------------------
# Optional: size distribution bins
# -----------------------------
if ("length" %in% names(mob_all)) {
  mob_all <- mob_all %>% mutate(length_num = suppressWarnings(as.numeric(length))) %>%
    mutate(size_bin = cut(length_num, breaks = c(-Inf,5e3,1e4,2e4,5e4,1e5,2e5,5e5,Inf),
                          labels = c('<5kb','5-10kb','10-20kb','20-50kb','50-100kb','100-200kb','200-500kb','>500kb')))
  size_dist <- mob_all %>% count(sample, size_bin)
  write_tsv(size_dist, "mobilome_size_bins.tsv")
}

# -----------------------------
# Session info for reproducibility
# -----------------------------
writeLines(c('==== SESSION INFO ====','', capture.output(sessionInfo())), con = "mobilome_session_info.txt")

message("Mobilome analysis complete.")
