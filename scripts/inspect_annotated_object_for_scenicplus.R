#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[[1]] == length(args)) default else args[[hit[[1]] + 1]]
}

project_dir <- Sys.getenv("PROJECT_DIR")
if (!nzchar(project_dir)) {
  stop("PROJECT_DIR is not set. Source project_env.sh before running this script.")
}
project_dir <- normalizePath(project_dir, mustWork = TRUE)
inputs_dir <- file.path(project_dir, "inputs")
results_dir <- file.path(project_dir, "results", "annotated_object")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- get_arg("--object-path", Sys.getenv("ANNOTATED_OBJECT"))
if (!nzchar(object_path)) {
  stop("Provide --object-path /path/to/object.rds or set ANNOTATED_OBJECT.")
}
object_path <- normalizePath(object_path, mustWork = TRUE)
object_format <- tolower(get_arg("--object-format", tools::file_ext(object_path)))

if (object_format %in% c("rds", "rda")) {
  obj <- readRDS(object_path)
} else if (object_format == "qs") {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("object_format=qs requires the R package qs in the SCENIC+ environment.")
  }
  obj <- qs::qread(object_path)
} else {
  stop("Supported annotated Seurat formats are rds and qs. Convert other objects to Seurat RDS/QS first.")
}
if (!inherits(obj, "Seurat")) {
  stop("The annotated object must be a Seurat object.")
}

meta <- obj@meta.data
cols <- colnames(meta)
score_col <- function(patterns) {
  score <- rep(0L, length(cols))
  names(score) <- cols
  lower <- tolower(cols)
  for (i in seq_along(patterns)) {
    score <- score + ifelse(grepl(patterns[[i]], lower), length(patterns) - i + 1L, 0L)
  }
  score
}
pick_col <- function(patterns, default_value = "") {
  score <- score_col(patterns)
  if (max(score) == 0) default_value else names(which.max(score))
}
sample_col <- pick_col(c("^sample$", "sample", "orig.ident", "library", "donor", "replicate"), "sample")
condition_col <- pick_col(c("^condition$", "condition", "genotype", "group", "treatment"), "condition")
requested_cell_label_column <- Sys.getenv("CELL_LABEL_COLUMN")
cell_label_column <- if (nzchar(requested_cell_label_column) && requested_cell_label_column %in% cols) {
  requested_cell_label_column
} else {
  pick_col(c("^cell_label$", "cell_type", "celltype", "lineage", "annotation", "cluster"), "cell_label")
}
barcode_col <- pick_col(c("^barcode$", "barcodes", "cell_barcode"), "")
original_cell_id_col <- pick_col(c("^original_cell_id$", "orig_cell", "cell_id"), "")
reduction_names <- names(obj@reductions)
pick_reduction <- function(names_vec) {
  if (length(names_vec) == 0) return("")
  lower <- tolower(names_vec)
  hit <- which(grepl("final.*wnn.*umap|wnn.*final.*umap|eyefinal.*wnn.*umap", lower))
  if (length(hit) > 0) return(names_vec[[hit[[1]]]])
  exact_priority <- c("wnn.umap", "rna.umap", "umap")
  for (candidate in exact_priority) {
    hit <- which(lower == candidate)
    if (length(hit) > 0) return(names_vec[[hit[[1]]]])
  }
  wnn_hits <- which(grepl("wnn.*umap|umap.*wnn", lower))
  if (length(wnn_hits) > 0) {
    active_hits <- wnn_hits[!grepl("whole|atlas|global|eyefinal|eyeinitial", lower[wnn_hits])]
    if (length(active_hits) > 0) return(names_vec[[active_hits[[1]]]])
    return(names_vec[[wnn_hits[[1]]]])
  }
  hit <- which(grepl("final.*rna.*umap|rna.*final.*umap", lower))
  if (length(hit) > 0) return(names_vec[[hit[[1]]]])
  hit <- which(grepl("rna.*umap|umap.*rna", lower))
  if (length(hit) > 0) return(names_vec[[hit[[1]]]])
  hit <- which(grepl("umap", lower))
  if (length(hit) > 0) return(names_vec[[hit[[1]]]])
  names_vec[[1]]
}
reduction <- pick_reduction(reduction_names)
assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
layers <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
layer <- if ("counts" %in% layers) "counts" else if (length(layers) > 0) layers[[1]] else "counts"

candidate_rows <- list(
  data.frame(section = "object", name = "n_cells", value = ncol(obj)),
  data.frame(section = "object", name = "n_features_default_assay", value = nrow(obj[[assay]])),
  data.frame(section = "assays", name = names(obj@assays), value = names(obj@assays)),
  data.frame(section = "reductions", name = names(obj@reductions), value = names(obj@reductions)),
  data.frame(section = "metadata_columns", name = cols, value = vapply(meta, function(x) {
    ux <- unique(as.character(x))
    paste(head(ux, 8), collapse = "|")
  }, character(1)))
)
candidates <- do.call(rbind, candidate_rows)
write.table(candidates, file.path(results_dir, "annotated_object_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

params <- data.frame(
  parameter = c(
    "object_path",
    "object_format",
    "assay",
    "layer",
    "sample_col",
    "condition_col",
    "cell_label_column",
    "barcode_col",
    "original_cell_id_col",
    "reduction"
  ),
  value = c(
    object_path,
    object_format,
    assay,
    layer,
    sample_col,
    condition_col,
    cell_label_column,
    barcode_col,
    original_cell_id_col,
    reduction
  ),
  stringsAsFactors = FALSE
)
write.table(params, file.path(inputs_dir, "annotated_object_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
selected <- data.frame(
  field = c("object_format", "assay", "layer", "sample_col", "condition_col", "cell_label_column", "barcode_col", "original_cell_id_col", "reduction"),
  selected_value = c(object_format, assay, layer, sample_col, condition_col, cell_label_column, barcode_col, original_cell_id_col, reduction),
  purpose = c(
    "Input object format detected from suffix or argument.",
    "RNA count assay to export.",
    "Raw-count layer to export.",
    "Biological sample identifier for matching fragments and pseudobulk grouping.",
    "Condition/group identifier used by downstream differential summaries.",
    "Metadata column exported as the standardized downstream cell_label.",
    "Raw barcode column matching the fragment-file barcode field, if present.",
    "Original Seurat cell-name column, if distinct from colnames(object).",
    "Embedding used for metacell construction and inspection."
  ),
  stringsAsFactors = FALSE
)
write.table(selected, file.path(results_dir, "annotated_object_selected_fields.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

label_preview <- if (cell_label_column %in% cols) {
  as.data.frame(sort(table(as.character(meta[[cell_label_column]])), decreasing = TRUE))
} else {
  data.frame()
}
if (nrow(label_preview) > 0) {
  colnames(label_preview) <- c("label", "n_cells")
  write.table(label_preview, file.path(results_dir, "annotated_object_label_preview.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
metadata_preview <- data.frame(
  column = cols,
  class = vapply(meta, function(x) paste(class(x), collapse = ","), character(1)),
  n_unique = vapply(meta, function(x) length(unique(as.character(x))), integer(1)),
  example_values = vapply(meta, function(x) {
    ux <- unique(as.character(x))
    paste(head(ux, 8), collapse = " | ")
  }, character(1)),
  stringsAsFactors = FALSE
)
write.table(metadata_preview, file.path(results_dir, "annotated_object_metadata_preview.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
report_path <- file.path(results_dir, "annotated_object_inspection_report.md")
report <- c(
  "# Annotated Object Inspection Report",
  "",
  "## Object",
  paste0("- Path: `", object_path, "`"),
  paste0("- Format: `", object_format, "`"),
  paste0("- Cells: ", ncol(obj)),
  paste0("- Selected assay: `", assay, "`"),
  paste0("- Selected layer: `", layer, "`"),
  paste0("- Assays: ", paste(names(obj@assays), collapse = ", ")),
  paste0("- Reductions: ", paste(names(obj@reductions), collapse = ", ")),
  "",
  "## Automatically Selected Fields",
  paste0("- sample_col: `", sample_col, "`"),
  paste0("- condition_col: `", condition_col, "`"),
  paste0("- CELL_LABEL_COLUMN: `", cell_label_column, "`"),
  paste0("- barcode_col: `", ifelse(nzchar(barcode_col), barcode_col, "<derive from cell names>"), "`"),
  paste0("- original_cell_id_col: `", ifelse(nzchar(original_cell_id_col), original_cell_id_col, "<object cell names>"), "`"),
  paste0("- reduction: `", reduction, "`"),
  "",
  "## Files Written",
  "- `inputs/annotated_object_params.tsv`",
  "- `results/annotated_object/annotated_object_selected_fields.tsv`",
  "- `results/annotated_object/annotated_object_metadata_preview.tsv`",
  "- `results/annotated_object/annotated_object_label_preview.tsv`"
)
writeLines(report, report_path)

cat("WROTE", file.path(inputs_dir, "annotated_object_params.tsv"), "\n")
cat("WROTE", file.path(results_dir, "annotated_object_candidates.tsv"), "\n")
cat("WROTE", file.path(results_dir, "annotated_object_selected_fields.tsv"), "\n")
cat("WROTE", file.path(results_dir, "annotated_object_metadata_preview.tsv"), "\n")
cat("WROTE", report_path, "\n")
if (nrow(label_preview) > 0) {
  cat("WROTE", file.path(results_dir, "annotated_object_label_preview.tsv"), "\n")
}
cat("\nANNOTATED OBJECT SUMMARY\n")
cat("cells:", ncol(obj), "\n")
cat("assay/layer:", assay, "/", layer, "\n")
cat("selected sample_col:", sample_col, "\n")
cat("selected condition_col:", condition_col, "\n")
cat("selected CELL_LABEL_COLUMN:", cell_label_column, "\n")
cat("selected reduction:", reduction, "\n")
cat("Review annotated_object_params.tsv before export. Confirm CELL_LABEL_COLUMN before continuing.\n")
