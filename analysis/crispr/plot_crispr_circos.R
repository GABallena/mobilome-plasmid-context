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

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CRISPR HOST↔PLASMID NETWORK VISUALIZATION - ULTRA EDITION                  ║
# ║  Next-gen circular network plot with cyberpunk aesthetics                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Update (2025-09-19): Added a publication-friendly 'scientific' theme as default,
# toned down effects, improved label safety, and added simple legends for clarity.

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ CONFIGURATION ZONE                                                          │
# └─────────────────────────────────────────────────────────────────────────────┘
set.seed(42)  # Reproducibility

paths <- list(
  spacer_hits    = "results/tables/spacer_edges_contig.tsv",
  contig_metrics = "results/tables/contig_metrics.tsv",
  plc_by_contig  = "results/tables/plc_class_by_contig.tsv",
  mob_root       = "mob_results"
)
 
# Visual parameters - high-res output
top_hosts    <- 50
top_plasmids <- 50
out_png      <- "results/figs/crispr_circos_ultra.png"
png_width    <- 4000
png_height   <- 4000
png_res      <- 300

# Package initialization
pkgs <- c("readr","dplyr","tidyr","stringr","purrr","circlize","optparse",
          "scales","ComplexHeatmap","viridis","colorspace","RColorBrewer")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(circlize); library(optparse); library(scales); library(ComplexHeatmap)
  library(grid); library(viridis); library(colorspace)

# ---- Config ----
EDGES_TSV <- Sys.getenv("EDGES_TSV", "results/tables/crispr_edges.tsv")
CONTIG_META_TSV <- Sys.getenv("CONTIG_META_TSV", "results/tables/contig_metadata.tsv")
OUT_DIR <- Sys.getenv("OUT_DIR", "results/figures/crispr")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

})

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ CITYCORE 12-Color Palette (only palette used)                              │
# └─────────────────────────────────────────────────────────────────────────────┘
# Background and base fills
bg_primary <- "white"
bg_glow <- scales::alpha("#000000", 0.02)

# CityCore 12 palette (from maps.R)
citystat_core12 <- c(
  Navy      = "#233B5D",
  Blue      = "#2F6DA8",
  DeepTeal  = "#1C7C7D",
  Cyan      = "#17BEBB",
  Forest    = "#2E7D32",
  Olive     = "#7A8F30",
  Mustard   = "#C7A41A",
  Orange    = "#E07A2D",
  Vermilion = "#D0432B",
  Crimson   = "#B01E2F",
  Plum      = "#7A3E9D",
  Indigo    = "#3B4BA3"
)

# Helpers to pick CityCore colors for semantic roles
cc <- list(
  host       = citystat_core12[["Forest"]],    # green
  plasmid    = citystat_core12[["Plum"]],      # purple
  accent     = citystat_core12[["Orange"]],    # orange
  gray       = "#636363",
  mobility_c = citystat_core12[["Forest"]],    # conjugative
  mobility_m = citystat_core12[["Blue"]],      # mobilizable
  mobility_n = "#636363",                       # non-mobilizable
  mobility_u = "#bdbdbd",                       # unknown
  relaxase   = citystat_core12[["Blue"]],
  orit       = citystat_core12[["Vermilion"]]
)

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ ENHANCED CLI INTERFACE                                                      │
# └─────────────────────────────────────────────────────────────────────────────┘
option_list <- list(
  optparse::make_option(c("-e","--edges"), type="character", default=paths$spacer_hits),
  optparse::make_option(c("--contigs"), type="character", default=paths$contig_metrics),
  optparse::make_option(c("--plc"), type="character", default=paths$plc_by_contig),
  optparse::make_option(c("--mob-dir"), type="character", default=paths$mob_root),
  optparse::make_option(c("--top-hosts"), type="integer", default=top_hosts),
  optparse::make_option(c("--top-plasmids"), type="integer", default=top_plasmids),
  optparse::make_option(c("--min-weight"), type="integer", default=1),
  # Threshold for drawing ribbons; edges below this are filtered out of the diagram
  optparse::make_option(c("--ribbon-threshold"), type="integer", default=1,
                       help="Minimum edge weight required to draw a ribbon"),
  optparse::make_option(c("--max-edges-per-host"), type="integer", default=8),
  optparse::make_option(c("-o","--out"), type="character", default=out_png),
  optparse::make_option(c("--theme"), type="character", default="scientific",
                       help="Visual theme: cyberpunk, vaporwave, aurora, matrix"),
  optparse::make_option(c("--glow-intensity"), type="double", default=0.8),
  optparse::make_option(c("--label-style"), type="character", default="holographic")
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))
theme_mode <- opt$theme
glow_intensity <- opt$`glow-intensity`

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ HELPER FUNCTIONS - ENHANCED                                                │
# └─────────────────────────────────────────────────────────────────────────────┘
base_contig <- function(x){
  x <- sub("\\s.*$", "", x)
  sub("_flag=.*$", "", x)
}

# Stylized labels with futuristic formatting
format_label <- function(id, style="cyber"){
  # For contigs without genus, extract just the k141 number
  if(grepl("k141_\\d+", id)) {
    # Extract just the contig number
    return(sub(".*k141_(\\d+).*", "Contig \\1", id))
  }
  
  # For plasmids without genus, show a cleaned version
  id <- sub("^PROJECT-", "", id)
  id <- sub("^[^-]*-", "", id)
  
  if(style == "cyber"){
    # Add subtle formatting
    parts <- strsplit(id, "_")[[1]]
    if(length(parts) > 2){
      id <- paste0("◊", parts[2], ".", tail(parts,1))
    }
  }
  id
}

# Glow effect simulator
add_glow <- function(color, intensity=0.5){
  scales::alpha(color, intensity * runif(1, 0.7, 1.0))
}

# Safe infix OR helper for defaults
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.atomic(a) && all(is.na(a))) return(b)
  a
}

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ DATA LOADING & PROCESSING                                                   │
# └─────────────────────────────────────────────────────────────────────────────┘
# Load PLSDB taxonomy mapping for plasmid genus names
plsdb_nuccore <- if(file.exists("PLSDB_db/nuccore.csv")) {
  readr::read_csv("PLSDB_db/nuccore.csv", show_col_types=FALSE)
} else {
  tibble(NUCCORE_ACC=character(), TAXONOMY_UID=character())
}

plsdb_taxonomy <- if(file.exists("PLSDB_db/taxonomy.csv")) {
  readr::read_csv("PLSDB_db/taxonomy.csv", show_col_types=FALSE)
} else {
  tibble(TAXONOMY_UID=character(), TAXONOMY_genus=character())
}

# Create accession to genus mapping for plasmids
acc2genus <- plsdb_nuccore %>%
  select(NUCCORE_ACC, TAXONOMY_UID) %>%
  left_join(plsdb_taxonomy %>% select(TAXONOMY_UID, TAXONOMY_genus), by="TAXONOMY_UID") %>%
  mutate(
    # Remove version numbers from accessions for matching
    acc_base = sub("\\.\\d+$", "", NUCCORE_ACC),
    genus_label = coalesce(TAXONOMY_genus, "Unknown")
  ) %>%
  select(acc_base, genus_label) %>%
  distinct()

# Load GTDB-Tk taxonomy for host contigs
gtdbtk_files <- list.files("gtdbtk_drep_out", pattern="gtdbtk.bac120.summary.tsv$", 
                           recursive=TRUE, full.names=TRUE)

if(length(gtdbtk_files) > 0) {
  gtdbtk_tax <- purrr::map_dfr(gtdbtk_files, ~{
    df <- readr::read_tsv(.x, show_col_types=FALSE)
    df %>% select(user_genome, classification)
  }) %>%
    distinct() %>%
    mutate(
      # Extract genus from GTDB classification string (g__GenusName)
      genus = str_match(classification, ";g__([^;]+)")[,2],
      genus = str_replace(genus, "_[A-Z]$", ""),  # Remove suffix like _A, _B
      # Filter out placeholder genome IDs (uppercase alphanumeric codes like JAJXWB01)
      genus = ifelse(str_detect(genus, "^[A-Z0-9]+$"), NA_character_, genus),
      genus = ifelse(is.na(genus) | genus == "", "Unknown", genus),
      mag_id = user_genome
    ) %>%
    select(mag_id, genus)
  
  cat(sprintf("Loaded %d MAG taxonomies from GTDB-Tk\n", nrow(gtdbtk_tax)))
  cat(sprintf("MAGs with proper genus names: %d\n", sum(gtdbtk_tax$genus != "Unknown")))
  cat("Sample MAG taxonomies:\n")
  print(head(gtdbtk_tax %>% filter(genus != "Unknown"), 5))
} else {
  gtdbtk_tax <- tibble(mag_id=character(), genus=character())
}

# Load contig-to-bin mapping
contig2bin <- if(file.exists("mag_analysis/full_contig_bin_map.tsv")) {
  readr::read_tsv("mag_analysis/full_contig_bin_map.tsv", show_col_types=FALSE) %>%
    mutate(
      # Match the exact format in spacer_edges_contig.tsv
      contig_full = paste0(sample, "_", contig)
    ) %>%
    select(contig_full, mag_id)
} else {
  tibble(contig_full=character(), mag_id=character())
}

# Create contig to genus mapping - map bin IDs to genus
contig2genus <- contig2bin %>%
  left_join(gtdbtk_tax, by="mag_id") %>%
  mutate(genus = coalesce(genus, "Unknown")) %>%
  select(contig_full, genus) %>%
  distinct()

# Debug: check if mapping worked
if(nrow(contig2genus) > 0) {
  cat(sprintf("Loaded %d contig-to-genus mappings\n", nrow(contig2genus)))
  cat("Sample mappings:\n")
  print(head(contig2genus, 3))
}

stopifnot(file.exists(paths$spacer_hits))
edges_raw <- readr::read_tsv(paths$spacer_hits, show_col_types=FALSE)

# Column standardization
nm <- names(edges_raw)
if("host_contig" %in% nm && !("contig" %in% nm)) {
  edges_raw <- dplyr::rename(edges_raw, contig=host_contig)
}
if("plasmid" %in% nm && !("plasmid_id" %in% nm)) {
  edges_raw <- dplyr::rename(edges_raw, plasmid_id=plasmid)
}

# Process edges
if("spacer_id" %in% names(edges_raw)) {
  edges <- edges_raw %>%
    mutate(contig=base_contig(contig)) %>%
    group_by(contig, plasmid_id) %>%
    summarise(weight=n_distinct(spacer_id), .groups="drop")
} else if("weight" %in% names(edges_raw)) {
  edges <- edges_raw %>%
    mutate(contig=base_contig(contig)) %>%
    group_by(contig, plasmid_id) %>%
    summarise(weight=sum(as.numeric(weight)), .groups="drop")
} else if("n_unique_spacers" %in% names(edges_raw)) {
  edges <- edges_raw %>%
    mutate(contig=base_contig(contig)) %>%
    group_by(contig, plasmid_id) %>%
    summarise(weight=sum(as.numeric(n_unique_spacers)), .groups="drop")
}

edges <- edges %>% filter(weight >= opt$`min-weight`)

# Load contig metrics
stopifnot(file.exists(paths$contig_metrics))
cm <- readr::read_tsv(paths$contig_metrics, show_col_types=FALSE) %>%
  mutate(length=as.numeric(length), gc=as.numeric(gc))

# MOB-suite annotations
cr_files <- if(dir.exists(paths$mob_root)) {
  list.files(paths$mob_root, pattern="contig_report.txt$", recursive=TRUE, full.names=TRUE)
} else character(0)

if(!length(cr_files)) {
  mob <- tibble(contig_id=character(), mobility=character(), has_relax=integer(), has_orit=integer())
} else {
  mob <- purrr::map_dfr(cr_files, ~{
    df <- suppressWarnings(readr::read_tsv(.x, show_col_types=FALSE))
    nm <- names(df)
    if("relaxase_type(s)" %in% nm) names(df)[match("relaxase_type(s)", nm)] <- "relaxase_types"
    if("orit_type(s)" %in% nm) names(df)[match("orit_type(s)", nm)] <- "orit_types"
    if(!("contig_id" %in% names(df))) {
      if("contig" %in% names(df)) names(df)[match("contig", names(df))] <- "contig_id"
      else if("id" %in% names(df)) names(df)[match("id", names(df))] <- "contig_id"
    }
    dplyr::select(df, dplyr::any_of(c("contig_id","predicted_mobility","relaxase_types","orit_types")))
  }) %>%
    mutate(contig_id=base_contig(contig_id),
           mobility=dplyr::coalesce(predicted_mobility, "-"),
           has_relax=as.integer(!is.na(relaxase_types) & relaxase_types!="-"),
           has_orit=as.integer(!is.na(orit_types) & orit_types!="-")) %>%
    select(contig_id, mobility, has_relax, has_orit) %>%
    distinct()
}

# PLC annotations
plc <- if(file.exists(paths$plc_by_contig)) {
  readr::read_tsv(paths$plc_by_contig, show_col_types=FALSE) %>%
    select(contig_id, plc_class)
} else tibble(contig_id=character(), plc_class=character())

# Merge annotations
host_annot <- edges %>%
  distinct(contig) %>%
  left_join(mob, by=c("contig"="contig_id")) %>%
  left_join(plc, by=c("contig"="contig_id")) %>%
  left_join(cm, by=c("contig"="contig_id")) %>%
  left_join(contig2genus, by=c("contig"="contig_full")) %>%
  mutate(
    mobility=dplyr::recode(mobility,
      "non-mobilizable"="non-mobilizable",
      "mobilizable"="mobilizable",
      "conjugative"="conjugative",
      .default="unknown"),
    gc=as.numeric(gc),
    length=as.numeric(length),
    genus = coalesce(genus, "Unknown")
  )

# Debug: check genus assignment
cat(sprintf("\nHost annotation summary:\n"))
cat(sprintf("Total hosts: %d\n", nrow(host_annot)))
cat(sprintf("Hosts with genus assigned: %d\n", sum(host_annot$genus != "Unknown")))
cat("Sample host annotations:\n")
print(host_annot %>% select(contig, genus) %>% head(5))

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ TOP NODE SELECTION WITH NETWORK METRICS                                    │
# └─────────────────────────────────────────────────────────────────────────────┘
host_rank <- edges %>% 
  group_by(contig) %>% 
  summarise(w=sum(weight), connections=n(), avg_weight=mean(weight), .groups="drop") %>% 
  arrange(desc(w))

plas_rank <- edges %>% 
  group_by(plasmid_id) %>% 
  summarise(w=sum(weight), connections=n(), avg_weight=mean(weight), .groups="drop") %>% 
  arrange(desc(w))

# Include ALL hosts that have a PLC assignment, and all their edges
hosts_with_plc <- host_annot %>% filter(!is.na(plc_class) & plc_class != "") %>% pull(contig)
host_order <- host_rank$contig[host_rank$contig %in% hosts_with_plc]

E <- edges %>% filter(contig %in% host_order)

# All plasmids connected to those hosts, ordered by total weight
plas_from_hosts <- unique(E$plasmid_id)
plas_order <- plas_rank$plasmid_id[plas_rank$plasmid_id %in% plas_from_hosts]
E <- E %>% filter(plasmid_id %in% plas_order)

# Do not cap edges per host; include all

host_annot <- host_annot %>% filter(contig %in% host_order)
plas_nodes <- tibble(plasmid_id=plas_order) %>%
  mutate(acc_base = sub("\\.\\d+$", "", plasmid_id)) %>%
  left_join(acc2genus, by="acc_base") %>%
  mutate(genus = coalesce(genus_label, "Unknown"))

# Apply ribbon threshold: remove edges with weight below threshold
if (!is.null(opt$`ribbon-threshold`) && is.finite(opt$`ribbon-threshold`)) {
  thr <- as.numeric(opt$`ribbon-threshold`)
  E <- E %>% filter(weight >= thr)
}

# Attach PLC class to edges (keep chromosomal; render transparently later)
host_plc_map <- host_annot %>% select(contig, plc_class)
E <- E %>% left_join(host_plc_map, by = "contig")

# Early exit if no ribbons remain
if (nrow(E) == 0) {
  dir.create(dirname(out_png), showWarnings=FALSE, recursive=TRUE)
  if(grepl("\\.pdf$", out_png, ignore.case=TRUE)) {
    cairo_pdf(out_png, width=png_width/72, height=png_height/72, bg=bg_primary)
  } else {
    png(out_png, width=png_width, height=png_height, res=png_res, bg=bg_primary)
  }
  par(bg=bg_primary, mar=c(2,2,2,2))
  plot.new()
  title(main=sprintf("CRISPR Host–Plasmid Network\n(no ribbons ≥ threshold, non-chromosomal)"), col.main="black")
  dev.off()
  quit(save="no")
}

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ ADVANCED COLOR MAPPING - NEON CYBERPUNK                                    │
# └─────────────────────────────────────────────────────────────────────────────┘

# Mobility colors
mob_cols <- c(
  conjugative      = cc$mobility_c,
  mobilizable      = cc$mobility_m,
  `non-mobilizable` = cc$mobility_n,
  unknown          = cc$mobility_u
)

# Plasmid sector color: use a single scientific color (no rainbow)
plasmid_color <- cc$plasmid

# GC content: plasma gradient
gc_col_fun <- circlize::colorRamp2(
  c(0.25, 0.40, 0.55, 0.70, 0.85),
  c(citystat_core12[["Navy"]], citystat_core12[["Indigo"]], citystat_core12[["Plum"]],
    citystat_core12[["Vermilion"]], citystat_core12[["Mustard"]])
)

# Length: deep ocean to surface gradient  
len_col_fun <- circlize::colorRamp2(
  log10(c(1e3, 1e4, 1e5, 1e6, 1e7)),
  c(citystat_core12[["DeepTeal"]], citystat_core12[["Blue"]], citystat_core12[["Cyan"]],
    citystat_core12[["Olive"]], citystat_core12[["Mustard"]])
)

## Plasmid connectivity coloring removed per request

# Sector colors with gradient effects
sectors_present <- unique(c(E$contig, E$plasmid_id))
host_order_present <- host_order[host_order %in% sectors_present]
plas_order_present <- plas_order[plas_order %in% sectors_present]
sector_order <- c(host_order_present, plas_order_present)

grid_col <- c(
  # Uniform host sector color (remove host mobility encoding)
  setNames(rep(cc$host, length(host_order_present)), host_order_present),
  # Plasmid sectors use single color
  setNames(rep(plasmid_color, length(plas_order_present)), plas_order_present)
)
grid_col <- grid_col[match(sector_order, names(grid_col))]

## Link colors: encode PLC class of the host (contig)
plc_colors <- list(
  chromosomal = citystat_core12[["Blue"]],
  putative_plasmid = citystat_core12[["Olive"]],
  mobilizable = citystat_core12[["Orange"]]
)

host_plc_map <- host_annot %>% select(contig, plc_class)
base_link_col <- vapply(E$plc_class, function(cls) {
  col <- plc_colors[[as.character(cls)]] %||% citystat_core12[["Olive"]]
  col
}, character(1))

# Weight-based alpha, but make chromosomal ribbons fully transparent
chr_mask <- !is.na(E$plc_class) & E$plc_class == "chromosomal"
alpha_vals <- scales::rescale(E$weight^0.6, to=c(0.35, 0.85))
alpha_vals[chr_mask] <- 0
E$link_col <- scales::alpha(base_link_col, alpha_vals)
E$link_border <- ifelse(chr_mask, NA, scales::alpha("#000000", 0.15))

# Dynamic line widths with weight emphasis
E$link_lwd <- scales::rescale(E$weight^0.7, to=c(0.6, 6))

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ VISUALIZATION ENGINE - ULTRA HIGH DEFINITION                                │
# └─────────────────────────────────────────────────────────────────────────────┘
dir.create(dirname(out_png), showWarnings=FALSE, recursive=TRUE)

# Initialize ultra-HD device
if(grepl("\\.pdf$", out_png, ignore.case=TRUE)) {
  cairo_pdf(out_png, width=png_width/72, height=png_height/72, bg=bg_primary)
} else {
  png(out_png, width=png_width, height=png_height, res=png_res, bg=bg_primary)
}

# Canvas setup with dark theme - expanded margins for outward labels
par(bg=bg_primary, col="black", col.axis="black", col.lab="black",
    col.main="black", fg="black", mar=c(8, 8, 6, 12))

# Circos parameters for smooth aesthetics
circos.clear()
# Dynamic gaps: keep total gaps << 360 to avoid allocation errors on large N
n_sec <- length(sector_order)
per_gap <- if (n_sec > 0) max(0.05, min(0.6, 360/(n_sec*8))) else 0.4
big_gap <- max(3, per_gap*10)
gap_after <- rep(per_gap, n_sec)
if(length(host_order_present) > 0 && length(plas_order_present) > 0){
  gap_after[length(host_order_present)] <- big_gap  # Host/plasmid separation
}

circos.par(
  track.margin=c(0.001, 0.001),
  gap.after=gap_after,
  cell.padding=c(0, 0, 0, 0),
  track.height=0.04,
  start.degree=90,  # Top start for dramatic effect
  clock.wise=TRUE
)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ MAIN CHORD DIAGRAM WITH ENHANCED VISUALS                                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
chordDiagram(
  x=E %>% select(contig, plasmid_id, value=weight),
  order=sector_order,
  grid.col=grid_col,
  directional=0,
  transparency=0.3,
  link.sort=TRUE,
  link.largest.ontop=TRUE,
  annotationTrack="grid",
  preAllocateTracks=4,
  col=E$link_col,
  link.lwd=pmax(1, E$link_lwd * 0.6),
  link.border=E$link_border
)

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ TRACK 1: SECTOR LABELS WITH GENUS NAMES                                    │
# └─────────────────────────────────────────────────────────────────────────────┘
circos.trackPlotRegion(track.index=1, ylim=c(0, 1), panel.fun=function(x, y) {
  si <- get.cell.meta.data("sector.index")
  xcenter <- get.cell.meta.data("xcenter")
  is_host <- si %in% host_annot$contig
  
  if(is_host) {
    # Get genus name for this host contig
    genus_name <- host_annot$genus[match(si, host_annot$contig)]
    # Only show label if we have a real genus name
    if(!is.na(genus_name) && genus_name != "Unknown") {
      lab <- genus_name
      col_lab <- "black"
      cex_lab <- 0.32
      circos.text(xcenter, 1.2, lab, facing="clockwise",
                  niceFacing=TRUE, cex=cex_lab, col=col_lab, font=2, adj=c(0, 0.5))
    }
  } else {
    # Plasmid labels: use genus name
    genus_name <- plas_nodes$genus[match(si, plas_nodes$plasmid_id)]
    # Only show label if we have a real genus name
    if(!is.na(genus_name) && genus_name != "Unknown") {
      lab <- genus_name
      col_lab <- "black"
      cex_lab <- 0.26
      circos.text(xcenter, 1.2, lab, facing="clockwise",
                  niceFacing=TRUE, cex=cex_lab, col=col_lab, adj=c(0, 0.5))
    }
  }
}, bg.border=NA)

## TRACK 2 and 3 removed (relaxase and oriT)

## TRACK 2: Mobility (outer)
circos.trackPlotRegion(track.index=2, panel.fun=function(x, y) {
  si <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  # Mask base so right side isn't grey when no host data
  circos.rect(xlim[1], 0.1, xlim[2], 0.9, col=bg_primary, border=NA)
  if(si %in% host_annot$contig){
    mob <- dplyr::coalesce(host_annot$mobility[match(si, host_annot$contig)], "unknown")
  col <- mob_cols[mob]
    circos.rect(xlim[1], 0.1, xlim[2], 0.9, col=scales::alpha(col, 0.85), border=NA)
  }
}, bg.border=NA)

## Tracks for GC (3) and Length (4) removed per request

## TRACK 5 removed (no longer exists); only 3 additional tracks now (2: mobility, 3: outer glow accounting label track is 1)

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ TRACK 6: PLC CLASS - CLASSIFICATION RING                                   │
# └─────────────────────────────────────────────────────────────────────────────┘
## TRACK 5 removed (PLC class ring). Ribbons will encode PLC class instead.

# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ TRACK 7: OUTER GLOW RING                                                   │
# └─────────────────────────────────────────────────────────────────────────────┘
circos.trackPlotRegion(track.index=3, ylim=c(0,1), bg.border=NA,
  panel.fun=function(x,y){
    si <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    
    # Determine intensity based on connections
    if(si %in% host_annot$contig){
      intensity <- host_rank$w[match(si, host_rank$contig)]
      intensity <- scales::rescale(log10(intensity+1), to=c(0.1, 0.8))
  col <- cc$host
  circos.rect(xlim[1], 0, xlim[2], 1, col=scales::alpha(col, 0.12), border=NA)
    } else if (si %in% plas_nodes$plasmid_id) {
  col <- cc$plasmid
  circos.rect(xlim[1], 0, xlim[2], 1, col=scales::alpha(col, 0.12), border=NA)
    } else {
  circos.rect(xlim[1], 0, xlim[2], 1, col=scales::alpha(bg_glow, 0.06), border=NA)
    }
  }
)

# Finalize plot: title and legends (publication-friendly)
try(title(main="CRISPR Host–Plasmid Network", col.main="black"), silent=TRUE)

# ---- LEGENDS ---------------------------------------------------------------
# Helpers for safe conditionals
.has <- function(x) !is.null(x) && length(x) > 0

## Removed host mobility legend per request
lg_mob <- NULL

## PLC legend is now for ribbon colors (not a track); only include classes present in E
plc_present <- intersect(unique(as.character(E$plc_class)), names(plc_colors))
plc_present <- setdiff(plc_present, "chromosomal")
if (length(plc_present) == 0) {
  lg_plc <- NULL
} else {
  lg_plc <- ComplexHeatmap::Legend(
    title = "Ribbons: PLC class",
    at = plc_present,
    legend_gp = grid::gpar(fill = unname(unlist(plc_colors[plc_present])), col = NA)
  )
}

## Removed GC and Length legends per request
lg_gc <- NULL
lg_len <- NULL

## Connectivity legend removed per request

# 3) Link encodings ----------------------------------------------------------
# Link color legend removed (neutral link color)
lg_link_col <- NULL

# Link width legend: map to spacer support (weight)
# We'll show 3 representative widths spanning min–max.
make_lwd_legend <- function(weights, lwds) {
  rng_w <- range(weights, na.rm = TRUE, finite = TRUE)
  if (!is.finite(diff(rng_w)) || diff(rng_w) == 0) {
    w_seq <- rep(ifelse(is.finite(rng_w[1]), rng_w[1], 1), 3)
  } else {
    w_seq <- pretty(rng_w, n = 3)
  }
  rng_lwd <- range(lwds, na.rm = TRUE, finite = TRUE)
  if (!all(is.finite(rng_lwd))) rng_lwd <- c(1, 6)
  lwd_seq <- scales::rescale(w_seq, to = rng_lwd)
  graphics_fns <- lapply(seq_along(w_seq), function(i) {
    lwd_i <- max(1, lwd_seq[i])
    function(x, y, w, h) {
      grid::grid.segments(x - 0.45*w, y, x + 0.45*w, y,
                          gp = grid::gpar(lwd = lwd_i, col = "black", lineend = "round"))
    }
  })
  ComplexHeatmap::Legend(
    title = "Link width\n(# unique spacers)",
    at = as.character(w_seq),
    type = "grob",
    graphics = graphics_fns,
    labels = as.character(w_seq)
  )
}
lg_link_w <- make_lwd_legend(E$weight, E$link_lwd)

## Molecular markers legend removed (tracks removed)
lg_markers <- NULL

# 5) Sector meaning -----------------------------------------------------------
# A tiny legend clarifying ring sides
lg_sectors <- ComplexHeatmap::Legend(
  title = "Sectors",
  at = c("Hosts (contigs)", "Plasmids"),
  legend_gp = grid::gpar(
    fill = c(scales::alpha(cc$host, 0.8),
             scales::alpha(cc$plasmid, 0.8)),
    col = NA
  )
)

## Pack left-side legends (PLC class for ribbons)
left_legends <- list(lg_plc)
left_legends <- Filter(.has, left_legends)
if (length(left_legends)) {
  leg_left <- do.call(ComplexHeatmap::packLegend,
                      c(left_legends, list(direction = "vertical", gap = unit(4, "mm"))))
  ComplexHeatmap::draw(leg_left, x = unit(0.015, "npc"), y = unit(0.5, "npc"),
                       just = c("left", "center"))
}

## Pack the remaining legends on the right
right_legends <- list(
  lg_sectors,
  lg_link_w,
  lg_markers
)
right_legends <- Filter(.has, right_legends)
if (length(right_legends)) {
  leg_right <- do.call(ComplexHeatmap::packLegend,
                       c(right_legends, list(direction = "vertical", gap = unit(4, "mm"))))
  ComplexHeatmap::draw(leg_right, x = unit(0.985, "npc"), y = unit(0.5, "npc"),
                       just = c("right", "center"))
}

# ---------------------------------------------------------------------------
try(circos.clear(), silent=TRUE)
try(dev.off(), silent=TRUE)

# Generate HD poster version at 900 DPI
message("\n=== Creating HD poster version (900 DPI) ===")
out_png_hd <- "HD_poster_figures/abstract6_CRISPR_circos_HD.png"
dir.create(dirname(out_png_hd), showWarnings = FALSE, recursive = TRUE)

png_width_hd <- 6000
png_height_hd <- 6000
png_res_hd <- 900

png(out_png_hd, width=png_width_hd, height=png_height_hd, res=png_res_hd, bg=bg_primary)
par(bg=bg_primary, mar=c(0,0,0,0))

# Re-run the same circos plot with HD settings
circos.clear()
circos.par(start.degree=90, gap.degree=1, cell.padding=c(0,0,0,0),
           track.margin=c(0.005,0.005), points.overflow.warning=FALSE)

circos.initialize(factors=all_ids, xlim=matrix(c(0, seg_len[all_ids]), ncol=2, byrow=TRUE))

circos.track(ylim=c(0,1), panel.fun=function(x,y){
  sector_i <- get.cell.meta.data("sector.index")
  xcenter <- get.cell.meta.data("xcenter")
  ycenter <- get.cell.meta.data("ycenter")
  is_host <- sector_i %in% host_ids
  is_plas <- sector_i %in% plas_ids
  
  col_fill <- if(is_host) host_colors[sector_i] else if(is_plas) plasmid_colors[sector_i] else bg_primary
  col_border <- if(is_host) host_border else if(is_plas) plasmid_border else accent_primary
  
  circos.rect(xleft=0, ybottom=0, xright=seg_len[sector_i], ytop=1, col=col_fill, border=col_border, lwd=lwd_arcs)
  
  if(sector_i %in% names(label_map)){
    lab_txt <- label_map[sector_i]
    if(is_host){
      circos.text(x=xcenter, y=0.5, labels=lab_txt, facing="bending.inside", niceFacing=TRUE,
                  col=text_primary, cex=label_cex_host, font=2)
    } else if(is_plas){
      circos.text(x=xcenter, y=0.5, labels=lab_txt, facing="bending.outside", niceFacing=TRUE,
                  col=text_secondary, cex=label_cex_plasmid, font=1)
    }
  }
}, bg.border=NA, track.height=track_ht)

for(i in seq_len(nrow(E))){
  h_id <- E$host_id[i]; p_id <- E$plasmid_id[i]
  if(!h_id %in% all_ids || !p_id %in% all_ids) next
  w_norm <- E$weight_norm[i]
  col_ribbon <- ribbon_colors[i]
  transparency_val <- 0.3 + 0.5*w_norm
  lwd_ribbon <- ribbon_lwd_min + (ribbon_lwd_max - ribbon_lwd_min)*w_norm
  circos.link(sector.index1=h_id, point1=c(0, seg_len[h_id]),
              sector.index2=p_id, point2=c(0, seg_len[p_id]),
              col=scales::alpha(col_ribbon, transparency_val), lwd=lwd_ribbon, border=NA)
}

title_y <- 1.12; title_cex <- 2.8
mtext(side=3, line=-2, text="CRISPR-Cas Spacer–Plasmid Network", cex=title_cex, col=text_primary, font=2, at=0.5, adj=0.5)
mtext(side=3, line=-4, text=paste0("Top ", top_hosts, " host genomes × Top ", top_plasmids, " plasmid contigs"),
      cex=1.2, col=text_secondary, font=1, at=0.5, adj=0.5)

try(circos.clear(), silent=TRUE)
dev.off()
message("Saved HD poster CRISPR circos: ", normalizePath(out_png_hd, winslash = "/", mustWork = FALSE))
