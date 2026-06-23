#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

project_dir <- Sys.getenv("PROJECT_DIR")
if (!nzchar(project_dir)) {
  stop("PROJECT_DIR is not set. Source project_env.sh before running this script.")
}
project_dir <- normalizePath(project_dir, mustWork = TRUE)
inputs_dir <- file.path(project_dir, "inputs")
work_dir <- file.path(project_dir, "work", "annotated_seurat")
results_dir <- file.path(project_dir, "results", "annotated_object")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

params_path <- file.path(inputs_dir, "annotated_object_params.tsv")
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
object_path <- param("object_path")
if (!nzchar(object_path)) {
  stop("annotated_object_params.tsv must define object_path")
}
if (!grepl("^/", object_path)) {
  object_path <- file.path(project_dir, object_path)
}
object_path <- normalizePath(object_path, mustWork = TRUE)
object_format <- tolower(param("object_format", tools::file_ext(object_path)))

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
sample_col <- param("sample_col", "sample")
condition_col <- param("condition_col", "condition")
cell_label_column <- param("cell_label_column", "")
if (!nzchar(cell_label_column)) {
  stop("annotated_object_params.tsv must define cell_label_column. Rerun inspect_annotated_object.py.")
}
barcode_col <- param("barcode_col", "")
original_cell_id_col <- param("original_cell_id_col", "")
assay <- param("assay", "RNA")
layer <- param("layer", "counts")
reduction <- param("reduction", "")

required_meta <- c(sample_col, condition_col, cell_label_column)
missing <- setdiff(required_meta, colnames(meta))
if (length(missing) > 0) {
  stop("Annotated object metadata lacks required columns: ", paste(missing, collapse = ", "))
}

cells_keep <- rownames(meta)
obj <- subset(obj, cells = cells_keep)
meta <- obj@meta.data

counts <- GetAssayData(obj, assay = assay, layer = layer)
if (ncol(counts) != nrow(meta) || !all(colnames(counts) == rownames(meta))) {
  stop("RNA count matrix columns do not match Seurat metadata rownames.")
}

sample_id <- as.character(meta[[sample_col]])
condition <- as.character(meta[[condition_col]])
cell_label <- as.character(meta[[cell_label_column]])
original_cell_id <- if (nzchar(original_cell_id_col) && original_cell_id_col %in% colnames(meta)) {
  as.character(meta[[original_cell_id_col]])
} else {
  colnames(counts)
}
barcode <- if (nzchar(barcode_col) && barcode_col %in% colnames(meta)) {
  as.character(meta[[barcode_col]])
} else {
  x <- original_cell_id
  for (sid in unique(sample_id)) {
    idx <- sample_id == sid
    x[idx] <- sub(paste0("^", sid, "_"), "", x[idx])
  }
  x
}
cell_id <- paste0(barcode, "-", sample_id)
if (anyDuplicated(cell_id) > 0) {
  stop("Duplicated cell_id values after barcode-sample_id conversion.")
}
colnames(counts) <- cell_id

cell_meta <- data.frame(
  cell_id = cell_id,
  original_cell_id = original_cell_id,
  barcode = barcode,
  sample_id = sample_id,
  condition = condition,
  cell_label = cell_label,
  source_label = cell_label,
  analysis_unit = "cell",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (col in c("genotype", "replicate", "nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "FRiP", "TSS.enrichment")) {
  if (col %in% colnames(meta)) {
    cell_meta[[col]] <- meta[[col]]
  }
}
if (nzchar(reduction) && reduction %in% names(obj@reductions)) {
  emb <- Embeddings(obj, reduction = reduction)
  emb <- emb[colnames(obj), , drop = FALSE]
  for (i in seq_len(ncol(emb))) {
    cell_meta[[paste0(reduction, "_", i)]] <- emb[, i]
  }
}

write.table(cell_meta, file.path(inputs_dir, "cell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = rownames(counts)), file.path(work_dir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(cell_id = colnames(counts)), file.path(work_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
Matrix::writeMM(counts, file.path(work_dir, "rna_counts.genes_by_cells.mtx"))
system2("gzip", c("-f", file.path(work_dir, "rna_counts.genes_by_cells.mtx")))

summary <- as.data.frame.matrix(table(cell_meta$sample_id, cell_meta$cell_label))
write.table(summary, file.path(results_dir, "cells_by_sample_and_label.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(
  data.frame(
    metric = c("cells", "genes", "labels", "samples"),
    value = c(ncol(counts), nrow(counts), length(unique(cell_label)), length(unique(sample_id)))
  ),
  file.path(results_dir, "annotated_object_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("WROTE", file.path(inputs_dir, "cell_metadata.tsv"), "\n")
cat("WROTE", file.path(work_dir, "rna_counts.genes_by_cells.mtx.gz"), "\n")
cat("cells", ncol(counts), "genes", nrow(counts), "\n")
