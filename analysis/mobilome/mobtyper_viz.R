#!/usr/bin/env Rscript

# mobtyper_viz.R — Visualization toolkit for MOB-suite MOB-typer outputs
#
# Features
# - Auto-discovers MOB-typer result files in a directory (multiple samples)
# - Robust parser for typical MOB-typer tabular outputs
# - Summaries and publication-ready plots for:
#     1) Replicon types (Inc/rep_cluster, etc.) — stacked bar + heatmap (Top-N)
#     2) Relaxase types (MOBP/MOBQ/...) — stacked bar
#     3) MPF types (e.g., MPF_T) — stacked bar
#     4) oriT types — stacked bar
#     5) Predicted mobility categories — bar chart
# - Exports CSV summaries and high-res PNG/PDF figures
#
# Usage
#   Rscript mobtyper_viz.R \
#     --input-dir mobtyper_results \
#     --output-dir mobtyper_viz \
#     --pattern "(.*)\\.(txt|tsv|csv)$" \
#     --sample-regex "^(.*?)(?:[._]mobtyper.*)?$" \
#     --top-n 40 \
#     --count-mode count \
#     --palette economist \
#     --heatmap-normalize none --heatmap-cluster TRUE
#
# Notes
# - Expects MOB-typer aggregated outputs with columns similar to:
#   sample_id, num_contigs, size, gc, md5, rep_type(s), rep_type_accession(s),
#   relaxase_type(s), relaxase_type_accession(s), mpf_type, mpf_type_accession(s),
#   orit_type(s), orit_accession(s), predicted_mobility, mash_nearest_neighbor,
#   mash_neighbor_distance, mash_neighbor_identification, primary_cluster_id,
#   secondary_cluster_id, predicted_host_range_overall_* ...

suppressWarnings({})
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(rlang)
  library(forcats)
  library(patchwork)
})

# Quiet R CMD check for NSE columns used in dplyr/ggplot2
utils::globalVariables(c("name", "tot", "name2", "n", "sample", "value"))

# -- Helpers -------------------------------------------------------------------

make_clean <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

first_present <- function(nms, candidates) {
  for (c in candidates) {
    if (c %in% nms) return(c)
  }
  NA_character_
}

infer_sample_name <- function(filepath, sample_regex) {
  base <- basename(filepath)
  base <- sub("\\.(txt|tsv|csv|gz)$", "", base, ignore.case = TRUE)
  default_pat <- "^(.*?)(?:[._]mobtyper.*)?$"
  pat <- if (is.character(sample_regex) && length(sample_regex) == 1L && !is.na(sample_regex) && nzchar(sample_regex)) sample_regex else default_pat
  m <- tryCatch(regexec(pat, base, perl = TRUE), error = function(e) regexec(default_pat, base, perl = TRUE))
  matched <- regmatches(base, m)
  if (length(matched) == 1 && length(matched[[1]]) >= 2) return(matched[[1]][2])
  sub("[._]mobtyper.*$", "", base, perl = TRUE)
}

save_plot <- function(p, file_base, width = 10, height = 6, dpi = 300, bg = "white",
                      png_width = NA_real_, png_height = NA_real_, png_dpi = NA_real_,
                      pdf_width = NA_real_, pdf_height = NA_real_) {
  # Allow per-call overrides first; if not provided, fall back to CLI overrides; else use plot width/height
  png_w_cli <- tryCatch(as.numeric(get("args", inherits = TRUE)$`png-width`), error = function(e) NA_real_)
  png_h_cli <- tryCatch(as.numeric(get("args", inherits = TRUE)$`png-height`), error = function(e) NA_real_)
  png_d_cli <- tryCatch(as.numeric(get("args", inherits = TRUE)$`png-dpi`), error = function(e) NA_real_)
  pdf_w_cli <- tryCatch(as.numeric(get("args", inherits = TRUE)$`pdf-width`), error = function(e) NA_real_)
  pdf_h_cli <- tryCatch(as.numeric(get("args", inherits = TRUE)$`pdf-height`), error = function(e) NA_real_)

  w_png <- if (!is.na(png_width)) png_width else if (!is.na(png_w_cli)) png_w_cli else width
  h_png <- if (!is.na(png_height)) png_height else if (!is.na(png_h_cli)) png_h_cli else height
  d_png <- if (!is.na(png_dpi)) png_dpi else if (!is.na(png_d_cli)) png_d_cli else dpi
  w_pdf <- if (!is.na(pdf_width)) pdf_width else if (!is.na(pdf_w_cli)) pdf_w_cli else width
  h_pdf <- if (!is.na(pdf_height)) pdf_height else if (!is.na(pdf_h_cli)) pdf_h_cli else height

  ggsave(paste0(file_base, ".png"), p, width = w_png, height = h_png, dpi = d_png, bg = bg, limitsize = FALSE)
  pdf_dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
  ggsave(paste0(file_base, ".pdf"), p, width = w_pdf, height = h_pdf, device = pdf_dev, bg = bg, limitsize = FALSE)
}

get_base_palette <- function(name = "economist") {
  name <- tolower(name)
  if (name == "ft" || name == "financial_times") {
    c("#173E6A", "#4899B8", "#B0133D", "#D65A76", "#E393A4", "#9DAA62", "#9D8D7B", "#C9B9A9")
  } else {
    c("#2F6FA1", "#41A9B8", "#3BB3DC", "#E2931B", "#C7B434", "#E5251F", "#E9541B", "#907070", "#C4B7A4")
  }
}

get_palette <- function(name = "economist", n = 10) {
  base <- get_base_palette(name)
  if (n <= length(base)) return(base[seq_len(n)])
  reps <- ceiling(n / length(base))
  cols <- rep(base, reps)[seq_len(n)]
  if (requireNamespace("colorspace", quietly = TRUE)) {
    cyc <- rep(seq_len(reps), each = length(base))[seq_len(n)]
    shift <- (cyc - 1) * 0.07
    cols <- colorspace::lighten(cols, amount = pmin(shift, 0.35))
  }
  cols
}

wrap_lab <- function(x, width = 18) {
  stringr::str_replace_all(stringr::str_wrap(x, width = width), "\n", "\n")
}

heatmap_scale_layer <- function(norm_mode) {
  low <- "#1B5F85"; mid <- "#F1EFEA"; high <- "#A61C2E"
  nm <- tolower(norm_mode)
  if (nm == "row_z") return(scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0, oob = scales::squish, name = "Row z-score"))
  scale_fill_gradient(low = low, high = high, oob = scales::squish, name = "Count")
}

normalize_matrix <- function(mat, mode = "none") {
  m <- as.matrix(mat)
  mode <- tolower(mode)
  if (mode == "row_z") {
    m <- t(scale(t(m)))
    m[is.nan(m)] <- 0
  } else if (mode == "row_prop") {
    rs <- rowSums(m); rs[rs == 0] <- 1; m <- m / rs
  } else if (mode == "col_prop") {
    cs <- colSums(m); cs[cs == 0] <- 1; m <- sweep(m, 2, cs, "/")
  } else if (mode == "log1p") {
    m <- log1p(m)
  }
  m
}

dist_matrix <- function(m, method = "euclidean") {
  method <- tolower(method)
  if (method %in% c("euclidean", "manhattan")) return(dist(m, method = method))
  if (method %in% c("pearson", "spearman")) {
    cmat <- suppressWarnings(stats::cor(t(m), method = method, use = "pairwise.complete.obs"))
    return(as.dist(1 - cmat))
  }
  dist(m)
}

cluster_order <- function(m, distance = "euclidean", linkage = "complete", axis = "rows") {
  if (axis == "cols") m <- t(m)
  d <- dist_matrix(m, distance)
  hc <- hclust(d, method = linkage)
  hc$order
}

# -- CLI -----------------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input-dir"), type = "character", default = "mobtyper_results", help = "Directory with MOB-typer files."),
  make_option(c("-o", "--output-dir"), type = "character", default = "mobtyper_viz", help = "Output directory for figures and CSVs."),
  make_option(c("-p", "--pattern"), type = "character", default = "(.*)\\.(txt|tsv|csv)$", help = "Regex for filenames to include."),
  make_option(c("-r", "--sample-regex"), type = "character", default = "^(.*?)(?:[._]mobtyper.*)?$", help = "Regex with one capture group to extract sample name from filename."),
  make_option(c("-n", "--top-n"), type = "integer", default = 40L, help = "Top N features for plots."),
  make_option(c("--count-mode"), type = "character", default = "count", help = "Counting mode for list columns: count|presence [default: %default]"),
  make_option(c("--palette"), type = "character", default = "economist", help = "Palette name: economist|ft [default: %default]"),
  make_option(c("--label-wrap"), type = "integer", default = 18L, help = "Label wrap width [default: %default]"),
  make_option(c("--heatmap-normalize"), type = "character", default = "none", help = "Heatmap normalization: none|row_z|row_prop|col_prop|log1p [default: %default]"),
  make_option(c("--heatmap-cluster"), type = "logical", default = TRUE, help = "Cluster rows/columns in heatmap [default: %default]"),
  make_option(c("--cluster-distance"), type = "character", default = "euclidean", help = "Distance for clustering: euclidean|manhattan|pearson|spearman [default: %default]"),
  make_option(c("--cluster-linkage"), type = "character", default = "complete", help = "Linkage for clustering: complete|average|ward.D2|single [default: %default]"),
  # Output sizing (defaults NA so other plots keep their intrinsic sizes; override per-plot as needed)
  make_option(c("--png-width"), type = "double", default = NA, help = "PNG width in inches (default: match plot)"),
  make_option(c("--png-height"), type = "double", default = NA, help = "PNG height in inches (default: match plot)"),
  make_option(c("--png-dpi"), type = "integer", default = NA, help = "PNG DPI (default: 300)"),
  make_option(c("--pdf-width"), type = "double", default = NA, help = "PDF width inches (default: match plot)"),
  make_option(c("--pdf-height"), type = "double", default = NA, help = "PDF height inches (default: match plot)"),
  make_option(c("--metadata-file"), type = "character", default = "PROJECT DATA MASTERLIST - PROJECT Y1 DATA.tsv", help = "Optional metadata file to infer site and treatment. [default: %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

dir.create(args$`output-dir`, showWarnings = FALSE, recursive = TRUE)

# -- IO discovery --------------------------------------------------------------

all_files <- list.files(args$`input-dir`, full.names = TRUE)
if (length(all_files) == 0) {
  stop("No files found in ", args$`input-dir`)
}
keep <- grepl(args$pattern, basename(all_files), perl = TRUE, ignore.case = TRUE)
files <- all_files[keep]
if (length(files) == 0) stop("No files matching pattern in ", args$`input-dir`)

# -- Parsing -------------------------------------------------------------------

read_mob_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  # Read all columns as character to avoid type conflicts across files
  ct <- readr::cols(.default = readr::col_character())
  df <- if (ext %in% c("tsv", "txt")) {
    suppressWarnings(readr::read_tsv(path, col_types = ct))
  } else {
    suppressWarnings(readr::read_csv(path, col_types = ct))
  }
  names(df) <- make_clean(names(df))
  df
}

split_list_col <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "-"]
  if (!length(x)) return(character(0))
  # Values are comma-separated; trim spaces
  toks <- unlist(strsplit(paste(x, collapse = ","), ",", perl = FALSE, fixed = TRUE))
  toks <- trimws(toks)
  toks <- toks[nzchar(toks) & toks != "-"]
  toks
}

count_list_col <- function(df, sample, col, mode = "count") {
  vals <- split_list_col(df[[col]])
  if (!length(vals)) return(tibble(name = character(0), n = integer(0), sample = character(0)))
  if (tolower(mode) == "presence") {
    vals <- unique(vals)
    tibble(name = vals, n = 1L, sample = sample)
  } else {
    as_tibble(table(vals), .name_repair = "minimal") %>%
      transmute(name = as.character(vals), n = as.integer(n)) %>%
      mutate(sample = sample)
  }
}

parse_one <- function(path) {
  df <- read_mob_file(path)
  sample <- infer_sample_name(path, args$`sample-regex`)
  df$source_file <- basename(path)
  df$sample_file <- sample
  # Prefer sample_file as canonical sample name
  df$sample <- df$sample_file
  df
}

raw_tbls <- lapply(files, parse_one)
mob <- bind_rows(raw_tbls)

# Canonical names of known columns
nm <- names(mob)
col_map <- list(
  sample_id = first_present(nm, c("sample_id")),
  rep_types = first_present(nm, c("rep_type_s", "rep_types")),
  relaxase = first_present(nm, c("relaxase_type_s", "relaxase_types")),
  mpf = first_present(nm, c("mpf_type", "mpf")),
  orit = first_present(nm, c("orit_type_s", "orit_types")),
  mobility = first_present(nm, c("predicted_mobility", "mobility"))
)

# -- Summaries -----------------------------------------------------------------

count_mode <- tolower(args$`count-mode`)

summarize_feature <- function(df, colname, sample, mode = count_mode) {
  if (is.na(colname)) return(tibble())
  count_list_col(df, sample, colname, mode)
}

by_sample <- mob %>% group_split(sample, .keep = TRUE)

rep_counts <- map_dfr(by_sample, function(d) summarize_feature(d, col_map$rep_types, unique(d$sample)[1], count_mode))
relax_counts <- map_dfr(by_sample, function(d) summarize_feature(d, col_map$relaxase, unique(d$sample)[1], count_mode))
mpf_counts <- map_dfr(by_sample, function(d) summarize_feature(d, col_map$mpf, unique(d$sample)[1], count_mode))
orit_counts <- map_dfr(by_sample, function(d) summarize_feature(d, col_map$orit, unique(d$sample)[1], count_mode))

mobility_counts <- tibble()
if (!is.na(col_map$mobility)) {
  mobility_counts <- mob %>%
    mutate(.mob = .data[[col_map$mobility]]) %>%
    mutate(.mob = ifelse(.mob %in% c("", "-"), NA, .mob)) %>%
    distinct(sample, .mob) %>%
    count(sample, .mob, name = "n") %>%
    rename(mobility = .mob)
}

# -- Write CSVs ---------------------------------------------------------------

outdir <- args$`output-dir`

write_wide_counts <- function(df, fname) {
  if (nrow(df) == 0) return(invisible(NULL))
  wide <- tidyr::pivot_wider(df, names_from = sample, values_from = n, values_fill = 0)
  readr::write_csv(wide, file.path(outdir, fname))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write_wide_counts(rep_counts, "replicon_counts_wide.csv")
write_wide_counts(relax_counts, "relaxase_counts_wide.csv")
write_wide_counts(mpf_counts, "mpf_counts_wide.csv")
write_wide_counts(orit_counts, "orit_counts_wide.csv")
if (nrow(mobility_counts)) readr::write_csv(mobility_counts, file.path(outdir, "mobility_counts.csv"))

# -- Plotting -----------------------------------------------------------------

pal <- get_palette(args$palette, 20)

stacked_bar <- function(df, title, top_n = args$`top-n`, fname, wrap = args$`label-wrap`) {
  if (nrow(df) == 0) return(invisible(NULL))
  # Top features across all samples
  top_feats <- df %>% group_by(name) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot)) %>% slice_head(n = top_n) %>% pull(name)
  plot_df <- df %>% mutate(name2 = ifelse(name %in% top_feats, name, "Other")) %>%
    group_by(sample, name2) %>% summarise(n = sum(n), .groups = "drop") %>%
    mutate(name2 = fct_reorder(name2, n, sum, .desc = TRUE))
  p <- ggplot(plot_df, aes(x = sample, y = n, fill = name2)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = get_palette(args$palette, nlevels(plot_df$name2))) +
    labs(x = "Sample", y = "Count", title = title, fill = "Type") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot(p, file.path(outdir, fname), width = 12, height = 6)
}

heatmap_counts <- function(df, title, top_n = args$`top-n`, fname, norm = args$`heatmap-normalize`, cluster = args$`heatmap-cluster`) {
  if (nrow(df) == 0) return(invisible(NULL))
  top_feats <- df %>% group_by(name) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot)) %>% slice_head(n = top_n) %>% pull(name)
  mat <- df %>% filter(name %in% top_feats) %>%
    tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$name; mat$name <- NULL
  m <- as.matrix(mat)
  if (nrow(m) == 0 || ncol(m) == 0) return(invisible(NULL))
  m_norm <- normalize_matrix(m, norm)
  r_ord <- seq_len(nrow(m_norm))
  c_ord <- seq_len(ncol(m_norm))
  if (isTRUE(cluster)) {
    r_ord <- cluster_order(m_norm, distance = args$`cluster-distance`, linkage = args$`cluster-linkage`, axis = "rows")
    c_ord <- cluster_order(m_norm, distance = args$`cluster-distance`, linkage = args$`cluster-linkage`, axis = "cols")
  }
  df_long <- as_tibble(m_norm[ r_ord, c_ord, drop = FALSE ], rownames = "name") %>%
    pivot_longer(-name, names_to = "sample", values_to = "value")
  p <- ggplot(df_long, aes(x = sample, y = name, fill = value)) +
    geom_tile() +
    heatmap_scale_layer(norm) +
    labs(x = "Sample", y = "Feature", title = title) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, file.path(outdir, fname), width = 12, height = 10)
}

# Replicon types
stacked_bar(rep_counts, "Replicon types (Top-N)", args$`top-n`, "replicon_stacked")
heatmap_counts(rep_counts, "Replicon types — Heatmap (Top-N)", args$`top-n`, "replicon_heatmap")

# Relaxase types
stacked_bar(relax_counts, "Relaxase types", args$`top-n`, "relaxase_stacked")

# MPF types
stacked_bar(mpf_counts, "MPF types", args$`top-n`, "mpf_stacked")

# oriT types
stacked_bar(orit_counts, "oriT types", args$`top-n`, "orit_stacked")

# Mobility categories
if (nrow(mobility_counts)) {
  p_mob <- mobility_counts %>%
    mutate(mobility = ifelse(is.na(mobility) | mobility == "", "Unknown", mobility)) %>%
    ggplot(aes(x = sample, y = n, fill = mobility)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = get_palette(args$palette, nlevels(factor(mobility_counts$mobility)))) +
    labs(x = "Sample", y = "Count", title = "Predicted mobility", fill = "Mobility") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_mob, file.path(outdir, "mobility_stacked"), width = 12, height = 6)
}

# -- Replicon types: richness line plot by treatment ---------------------------

# Compute presence (unique types) per sample from rep_counts
rep_presence <- rep_counts %>%
  group_by(sample, name) %>% summarise(present = as.integer(sum(n) > 0), .groups = "drop")
rep_richness <- rep_presence %>% group_by(sample) %>% summarise(richness = sum(present), .groups = "drop")

# Attach sample_code (e.g., SAMPLE-01 from SAMPLE-01_S1)
rep_richness <- rep_richness %>% mutate(sample_code = sub("_S.*$", "", sample))

# Try to read metadata to get site and treatment
meta_path <- args$`metadata-file`
site_map <- NULL
if (!is.null(meta_path) && nzchar(meta_path) && file.exists(meta_path)) {
  meta <- suppressWarnings(readr::read_tsv(meta_path, show_col_types = FALSE, progress = FALSE))
  names(meta) <- make_clean(names(meta))
  sc_col <- first_present(names(meta), c("sample_code", "sample", "code"))
  sd_col <- first_present(names(meta), c("sample_description", "description", "sample_desc"))
  st_col <- first_present(names(meta), c("sample_type", "type"))
  if (!is.na(sc_col)) {
    site_map <- meta %>% transmute(
      sample_code = .data[[sc_col]],
      sample_description = if (!is.na(sd_col)) .data[[sd_col]] else NA_character_,
      sample_type = if (!is.na(st_col)) .data[[st_col]] else NA_character_
    ) %>%
      mutate(
        # Determine treatment from type or description
        treatment = case_when(
          # Broader matching; check 'untreat' before 'treat'
          grepl("(?i)untreat", sample_type, perl = TRUE) ~ "untreated",
          grepl("(?i)treat", sample_type, perl = TRUE) ~ "treated",
          grepl("(?i)untreat", sample_description, perl = TRUE) ~ "untreated",
          grepl("(?i)treat", sample_description, perl = TRUE) ~ "treated",
          TRUE ~ NA_character_
        ),
        # Exclude controls/optimization
        is_bad = grepl("control|optimization", paste(sample_type, sample_description), ignore.case = TRUE),
        # Extract site name by stripping trailing ' treated/untreated [rep]'
        site = ifelse(!is.na(sample_description) & nzchar(sample_description),
                      trimws(gsub("(?i)\\s+(treated|untreated).*$", "", sample_description, perl = TRUE)),
                      sample_code),
        # Extract replicate # if present (e.g., Untreated 2)
        replicate_num = {
          mm <- tryCatch(stringr::str_match(sample_description, "(?i)(?:un)?treated\\D*(\\d+)"), error = function(e) matrix(NA_character_, nrow = length(sample_description), ncol = 2))
          suppressWarnings(as.integer(mm[,2]))
        }
      ) %>%
      filter(!is_bad) %>%
      select(sample_code, site, treatment, replicate_num)
  }
}

# If no metadata, infer treatment/site minimally from sample_code
if (is.null(site_map)) {
  site_map <- rep_richness %>% transmute(sample_code, site = sample_code, treatment = NA_character_, replicate_num = NA_integer_)
}

rr <- rep_richness %>% left_join(site_map, by = "sample_code")

# Select replicate 2 when there are duplicates per site+treatment
rr_sel <- rr %>%
  group_by(site, treatment) %>%
  mutate(rep_rank = ifelse(!is.na(replicate_num), replicate_num, 1L)) %>%
  arrange(desc(rep_rank), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# Prepare for line plot: two lines (treated vs untreated) across sites
rr_long <- rr_sel %>% filter(!is.na(site) & !is.na(treatment)) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "treated")))

if (nrow(rr_long)) {
  # Export CSV for transparency/debugging
  readr::write_csv(rr_long %>% select(site, treatment, sample, richness, replicate_num), file.path(outdir, "replicon_richness_by_site.csv"))

  # Keep only sites with both treatments for dumbbell
  both_sites <- rr_long %>% distinct(site, treatment) %>% count(site) %>% filter(n >= 2) %>% pull(site)
  rr_db <- rr_long %>% filter(site %in% both_sites)
  if (nrow(rr_db) > 0) {
    # Order by mean richness (ascending)
    site_order <- rr_db %>% group_by(site) %>% summarise(avg = mean(richness), .groups = "drop") %>% arrange(avg) %>% pull(site)

    # Build anonymized labels from sample IDs (sample_code) per site
    lab_df <- rr_db %>% group_by(site) %>% summarise(
      untreated_code = dplyr::first(sample_code[treatment == "untreated"]),
      treated_code   = dplyr::first(sample_code[treatment == "treated"]),
      .groups = "drop"
    ) %>% mutate(
      site_lab = dplyr::case_when(
        !is.na(untreated_code) & !is.na(treated_code) ~ paste0(untreated_code, " | ", treated_code),
        !is.na(untreated_code) ~ untreated_code,
        !is.na(treated_code) ~ treated_code,
        TRUE ~ NA_character_
      )
    )
    site_lab_levels <- lab_df %>% filter(site %in% site_order) %>% arrange(match(site, site_order)) %>% pull(site_lab)

    # Compute gap per site (treated - untreated)
    gap_df <- rr_db %>% select(site, treatment, richness) %>%
      tidyr::pivot_wider(names_from = treatment, values_from = richness) %>%
      left_join(lab_df %>% select(site, site_lab), by = "site") %>%
      mutate(gap = treated - untreated,
             max_val = pmax(treated, untreated, na.rm = TRUE),
             gap_party_max = ifelse(treated >= untreated, "T", "U"))

    # Join max/gap and labels back to long
    rr_db2 <- rr_db %>% left_join(gap_df %>% select(site, gap, max_val, site_lab), by = "site") %>%
      mutate(is_max = richness == max_val,
             site_lab = factor(site_lab, levels = site_lab_levels))

  # Nudge settings (template-like spacing)
  nudge_value <- 2

    # Main dumbbell with callouts and inline legend labels
    p_main <- rr_db2 %>%
      mutate(
        label_x = richness + if_else(is_max, nudge_value, -nudge_value),
        label_hjust = if_else(is_max, 0, 1)
      ) %>%
      ggplot(aes(x = richness, y = site_lab)) +
  geom_line(aes(group = site), color = "#E7E7E7", linewidth = 3.5) +
  geom_point(aes(color = treatment), size = 3) +
      # value labels
  geom_text(aes(x = label_x, label = richness, color = treatment, hjust = label_hjust),
        size = 3.25, show.legend = FALSE) +
      # text legend on row with largest gap
      geom_text(data = rr_db2 %>% group_by(site) %>% summarise(gap = first(gap), .groups = "drop") %>% top_n(1, wt = abs(gap)) %>%
                  left_join(rr_db2, by = "site"),
                aes(label = treatment, color = treatment),
                nudge_y = 0.5, fontface = "bold", size = 3.25, show.legend = FALSE) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "#989898"),
        axis.title = element_blank(),
        panel.grid = element_blank()
      ) +
  scale_color_manual(values = c(untreated = "#BF2F24", treated = "#2CA02C")) +
  scale_y_discrete(limits = site_lab_levels, expand = expansion(mult = c(0,0), add = c(0.3, 1.2)))

    # Gap panel
  df_gap <- gap_df %>% mutate(
      gap_label = paste0(ifelse(gap >= 0, "+", ""), abs(gap), gap_party_max)
  ) %>% mutate(site_lab = factor(site_lab, levels = site_lab_levels))

    p_gap <- df_gap %>%
      ggplot(aes(x = gap, y = site_lab)) +
  geom_text(aes(x = 0, label = gap_label, color = gap_party_max),
        fontface = "bold", size = 3.25) +
      # header "Diff" above the top row without altering limits
      annotate("text", x = 0, y = tail(site_lab_levels, 1), label = "Diff",
               vjust = -1.2, fontface = "bold", size = 3.25) +
      theme_void() +
      coord_cartesian(xlim = c(-0.05, 0.05)) +
      theme(
        plot.margin = margin(l = 0, r = 0, b = 0, t = 0),
        panel.background = element_rect(fill = "#EFEFE3", color = "#EFEFE3"),
        legend.position = "none"
      ) +
      scale_color_manual(values = c(T = "#2CA02C", U = "#BF2F24")) +
      scale_y_discrete(limits = site_lab_levels, expand = expansion(mult = c(0,0), add = c(0.3, 1.2)))

    # Combine using patchwork (two-column layout)
  p_whole <- p_main + p_gap + plot_layout(widths = c(1, 0.12))

    # Save
  # Only shrink the Cleveland plot PNG by default; PDFs keep full size
  save_plot(p_whole, file.path(outdir, "replicon_types_cleveland"), width = 12, height = 8,
            png_width = ifelse(is.na(args$`png-width`), 6, args$`png-width`),
            png_height = ifelse(is.na(args$`png-height`), 3, args$`png-height`),
            png_dpi = ifelse(is.na(args$`png-dpi`), 300, args$`png-dpi`))
  }
}

message("mobtyper_viz: Done. Outputs in ", normalizePath(outdir, mustWork = FALSE))
