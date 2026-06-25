#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat(
    "Usage: prepare_metacell_inputs_from_seurat.R\n\n",
    "Create metacells from an annotated Seurat object using parameters in\n",
    "$PROJECT_DIR/inputs/metacell_params.tsv. Writes metacell metadata,\n",
    "membership, RNA counts and summary files for downstream SCENIC+.\n",
    sep = ""
  )
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(Matrix)
})

project_dir <- Sys.getenv("PROJECT_DIR")
if (!nzchar(project_dir)) {
  stop("PROJECT_DIR is not set. Source project_env.sh before running this script.")
}
project_dir <- normalizePath(project_dir, mustWork = TRUE)
inputs_dir <- file.path(project_dir, "inputs")
work_dir <- file.path(project_dir, "work", "metacells")
results_dir <- file.path(project_dir, "results", "metacells")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

write_project_setting <- function(key, value) {
  config <- file.path(project_dir, "scenicplus_project.env")
  lines <- if (file.exists(config)) readLines(config, warn = FALSE) else character(0)
  assignment <- paste0(key, "='", gsub("'", "'\"'\"'", value), "'")
  idx <- grep(paste0("^", key, "="), lines)
  if (length(idx) > 0) {
    lines[[idx[[1]]]] <- assignment
  } else {
    lines <- c(lines, assignment)
  }
  writeLines(lines, config)
}

params_path <- file.path(inputs_dir, "metacell_params.tsv")
if (!file.exists(params_path)) {
  stop("Missing ", params_path)
}
params_df <- read.delim(params_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("parameter", "value") %in% colnames(params_df))) {
  stop(params_path, " must contain columns: parameter and value")
}
params <- setNames(params_df$value, params_df$parameter)
param <- function(key, default = "") {
  value <- params[[key]]
  if (is.null(value) || is.na(value) || !nzchar(value)) default else value
}
param_int <- function(key, default) as.integer(param(key, as.character(default)))
param_vec <- function(key, default) {
  value <- param(key, default)
  as.integer(trimws(strsplit(value, "[, ]+")[[1]]))
}

seurat_rds <- param("seurat_rds")
if (!nzchar(seurat_rds)) {
  stop("metacell_params.tsv must define seurat_rds")
}
if (!grepl("^/", seurat_rds)) {
  seurat_rds <- file.path(project_dir, seurat_rds)
}
seurat_rds <- normalizePath(seurat_rds, mustWork = TRUE)

cell_meta_path <- file.path(inputs_dir, "cell_metadata.tsv")
if (!file.exists(cell_meta_path)) {
  stop("Missing ", cell_meta_path, ". Create single-cell metadata before metacell aggregation.")
}
single_cell_backup <- file.path(inputs_dir, "cell_metadata.single_cell.tsv")
cell_meta_input_path <- if (file.exists(single_cell_backup)) single_cell_backup else cell_meta_path
cell_meta <- read.delim(cell_meta_input_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!file.exists(single_cell_backup) &&
    "analysis_unit" %in% colnames(cell_meta) &&
    any(as.character(cell_meta$analysis_unit) == "metacell")) {
  stop(
    "inputs/cell_metadata.tsv appears to contain metacells, but ",
    single_cell_backup,
    " is missing. Rerun Step 2 to restore single-cell metadata before rebuilding metacells."
  )
}
required_cols <- c("cell_id", "barcode", "sample_id", "cell_label")
missing_cols <- setdiff(required_cols, colnames(cell_meta))
if (length(missing_cols) > 0) {
  stop("cell_metadata.tsv missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (!"original_cell_id" %in% colnames(cell_meta)) {
  cell_meta$original_cell_id <- paste0(cell_meta$sample_id, "_", cell_meta$barcode)
}
if (!"condition" %in% colnames(cell_meta)) {
  cell_meta$condition <- "condition_1"
}

object_format <- tolower(tools::file_ext(seurat_rds))
if (object_format == "rds") {
  obj <- readRDS(seurat_rds)
} else if (object_format == "qs") {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("metacell_params.tsv points to a .qs object, but the R package qs is not installed. Convert the object to .rds or install qs in the SCENIC+ environment.")
  }
  obj <- qs::qread(seurat_rds)
} else {
  stop("metacell_params.tsv seurat_rds must be a Seurat .rds or .qs file.")
}
if (!all(cell_meta$original_cell_id %in% colnames(obj))) {
  missing <- setdiff(cell_meta$original_cell_id, colnames(obj))
  stop("Seurat object lacks ", length(missing), " cells listed in cell_metadata.tsv. Example: ", missing[[1]])
}
obj <- subset(obj, cells = cell_meta$original_cell_id)
cell_meta <- cell_meta[match(colnames(obj), cell_meta$original_cell_id), , drop = FALSE]

obj$scenicplus_cell_id <- cell_meta$cell_id
obj$scenicplus_sample_id <- cell_meta$sample_id
obj$scenicplus_cell_label <- cell_meta$cell_label
obj$scenicplus_condition <- cell_meta$condition

wgcna_name <- param("wgcna_name", "scenicplus_metacells")
assay <- param("assay", "RNA")
layer <- param("layer", "counts")
reduction <- param("reduction", "wnn.umap")
dims <- param_vec("dims", "1 2")
k <- param_int("k", 12)
max_shared <- param_int("max_shared", 2)
min_cells <- param_int("min_cells", 1)
target_metacells <- param_int("target_metacells", 1000000)
mode <- param("mode", "average")

if (!reduction %in% names(obj@reductions)) {
  stop("Reduction not found in Seurat object: ", reduction)
}
obj <- SetupForWGCNA(obj, wgcna_name = wgcna_name)
obj <- MetacellsByGroups(
  obj,
  group.by = c("scenicplus_sample_id", "scenicplus_cell_label"),
  ident.group = "scenicplus_cell_label",
  k = k,
  reduction = reduction,
  dims = dims,
  assay = assay,
  layer = layer,
  mode = mode,
  cells.use = colnames(obj),
  min_cells = min_cells,
  max_shared = max_shared,
  target_metacells = target_metacells,
  verbose = TRUE,
  wgcna_name = wgcna_name
)

mc <- GetMetacellObject(obj, wgcna_name = wgcna_name)
counts <- GetAssayData(mc, assay = assay, layer = layer)
old_mc_ids <- colnames(counts)
metacell_barcode <- sprintf("MC%06d", seq_along(old_mc_ids))
metacell_id <- paste0(metacell_barcode, "-", as.character(mc$scenicplus_sample_id))
if (anyDuplicated(metacell_id) > 0) {
  stop("Duplicated metacell_id values after naming.")
}
colnames(counts) <- metacell_id

source_meta <- cell_meta
source_by_original <- split(source_meta, source_meta$original_cell_id)
embedding_cols <- names(source_meta)[startsWith(names(source_meta), paste0(reduction, "_"))]
membership_rows <- list()
metacell_rows <- list()
for (i in seq_along(old_mc_ids)) {
  old_id <- old_mc_ids[[i]]
  merged <- as.character(mc@meta.data[old_id, "cells_merged"])
  originals <- trimws(strsplit(merged, ",")[[1]])
  originals <- originals[nzchar(originals)]
  matched <- do.call(rbind, source_by_original[originals])
  if (is.null(matched) || nrow(matched) == 0) {
    stop("No source cells matched for metacell ", old_id)
  }
  sample_values <- unique(as.character(matched$sample_id))
  label_values <- unique(as.character(matched$cell_label))
  condition_values <- unique(as.character(matched$condition))
  if (length(sample_values) != 1 || length(label_values) != 1) {
    stop("Metacell crosses sample or cell_label boundaries: ", old_id)
  }
  condition_value <- if (length(condition_values) == 1) condition_values[[1]] else paste(condition_values, collapse = ";")
  membership_rows[[i]] <- data.frame(
    metacell_id = metacell_id[[i]],
    metacell_barcode = metacell_barcode[[i]],
    original_metacell_id = old_id,
    original_cell_id = matched$original_cell_id,
    cell_id = matched$cell_id,
    barcode = matched$barcode,
    sample_id = matched$sample_id,
    cell_label = matched$cell_label,
    condition = matched$condition,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metacell_row <- data.frame(
    cell_id = metacell_id[[i]],
    barcode = metacell_barcode[[i]],
    sample_id = sample_values[[1]],
    condition = condition_value,
    cell_label = label_values[[1]],
    original_metacell_id = old_id,
    n_member_cells = nrow(matched),
    analysis_unit = "metacell",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (length(embedding_cols) >= 2) {
    emb <- matched[, embedding_cols, drop = FALSE]
    emb <- as.data.frame(lapply(emb, function(x) as.numeric(as.character(x))), check.names = FALSE)
    metacell_row <- cbind(metacell_row, as.data.frame(as.list(colMeans(emb, na.rm = TRUE)), check.names = FALSE))
  }
  metacell_rows[[i]] <- metacell_row
}
membership <- do.call(rbind, membership_rows)
metacell_meta <- do.call(rbind, metacell_rows)

if (!file.exists(single_cell_backup)) {
  file.copy(cell_meta_input_path, single_cell_backup)
}
write.table(membership, file.path(inputs_dir, "metacell_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metacell_meta, cell_meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metacell_meta, file.path(inputs_dir, "metacell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(gene = rownames(counts)), file.path(work_dir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(cell_id = colnames(counts)), file.path(work_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
Matrix::writeMM(counts, file.path(work_dir, "rna_counts.genes_by_metacells.mtx"))
system2("gzip", c("-f", file.path(work_dir, "rna_counts.genes_by_metacells.mtx")))

summary <- as.data.frame.matrix(table(metacell_meta$sample_id, metacell_meta$cell_label))
write.table(summary, file.path(results_dir, "metacells_by_sample_and_label.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(
  data.frame(
    metric = c(
      "source_cells_in_metadata",
      "source_cells_in_metacells_unique",
      "membership_rows",
      "metacells",
      "median_member_cells",
      "max_member_cells",
      "max_shared_allowed",
      "k"
    ),
    value = c(
      nrow(cell_meta),
      length(unique(membership$cell_id)),
      nrow(membership),
      nrow(metacell_meta),
      median(metacell_meta$n_member_cells),
      max(metacell_meta$n_member_cells),
      max_shared,
      k
    )
  ),
  file.path(results_dir, "metacell_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("WROTE", file.path(inputs_dir, "metacell_membership.tsv"), "\n")
cat("WROTE", cell_meta_path, "\n")
cat("WROTE", file.path(work_dir, "rna_counts.genes_by_metacells.mtx.gz"), "\n")
cat("metacells", nrow(metacell_meta), "source_cells", nrow(cell_meta), "\n")
write_project_setting("ACTIVE_CELL_METADATA", cell_meta_path)
write_project_setting("METACELL_MEMBERSHIP", file.path(inputs_dir, "metacell_membership.tsv"))
write_project_setting("METACELL_SUMMARY", file.path(results_dir, "metacell_summary.tsv"))
