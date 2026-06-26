#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  if (any(x %in% c("--help", "-h"))) {
    cat("Usage: spgrn-plot-scenicplus-condition-stats [--outdir DIR] [--file-suffix SUFFIX] [--plot-style-config TSV] [--priority-eregulons TSV]\n")
    quit(status = 0)
  }
  out <- list(
    outdir = "results/scenicplus_figures",
    file_suffix = "",
    plot_style_config = "results/scenicplus_figures/plot_style_parameters.tsv",
    priority_eregulons = ""
  )
  i <- 1
  while (i <= length(x)) {
    key <- gsub("-", "_", sub("^--", "", x[[i]]))
    if (i == length(x) || grepl("^--", x[[i + 1]])) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- x[[i + 1]]
      i <- i + 2
    }
  }
  out
}

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), mustWork = FALSE)
project_path <- function(x) if (grepl("^/", x)) normalizePath(x, mustWork = FALSE) else normalizePath(file.path(project_dir, x), mustWork = FALSE)
named <- function(stem, suffix, ext = ".pdf") if (suffix == "") paste0(stem, ext) else paste0(stem, "_", suffix, ext)

args <- parse_args()
outdir <- project_path(args$outdir)
suffix <- args$file_suffix
style_path <- project_path(args$plot_style_config)
style_defaults <- data.frame(
  parameter = c(
    "font_family", "base_size", "axis_text_size", "axis_title_size", "plot_title_size", "strip_text_size",
    "legend_text_size", "legend_title_size", "caption_size", "axis_line_width", "heatmap_tile_border_width",
    "dot_min_size", "dot_size_scale", "umap_point_size", "umap_point_alpha", "embedding_point_size",
    "embedding_point_alpha", "condition_background_alpha", "condition_background_point_size", "max_panels_per_row",
    "color_control", "color_comparison",
    "color_group_1", "color_group_2", "color_group_3", "color_group_4", "color_group_5", "color_group_6"
  ),
  value = c(
    "Arial", "6.5", "5.6", "6.3", "7", "6.2",
    "5.5", "5.8", "5.0", "0.30", "0",
    "0.7", "3.0", "0.20", "0.75", "0.25",
    "0.75", "0.18", "0.10", "3",
    "#2F4858", "#B57457",
    "#2F4858", "#B57457", "#667A9A", "#9A8067", "#7D9477", "#8E5A63"
  ),
  stringsAsFactors = FALSE
)
style_info <- data.frame(
  parameter = style_defaults$parameter,
  meaning = c(
    "Font family used by Cairo PDF when available.",
    "Base ggplot font size in points.",
    "Axis tick label font size.",
    "Axis title font size.",
    "Panel title font size.",
    "Facet strip label font size.",
    "Legend item font size.",
    "Legend title font size.",
    "Figure legend/caption font size.",
    "Axis and tick line width.",
    "Heatmap tile border width; keep 0 for dense heatmaps.",
    "Minimum dot size for dot heatmaps.",
    "Maximum dot size for dot heatmaps.",
    "Point size for eRegulon AUC embedding overlays.",
    "Point alpha for eRegulon AUC embedding overlays.",
    "Point size for eRegulon activity UMAP/embedding panels.",
    "Point alpha for eRegulon activity UMAP/embedding panels.",
    "Background cell alpha in condition-resolved embedding panels.",
    "Background cell point size in condition-resolved embedding panels.",
    "Maximum number of condition/facet panels per row.",
    "Reference/control-like condition color when detected.",
    "Comparison/non-reference condition color.",
    "Generic categorical color 1 for cell labels or groups.",
    "Generic categorical color 2 for cell labels or groups.",
    "Generic categorical color 3 for cell labels or groups.",
    "Generic categorical color 4 for cell labels or groups.",
    "Generic categorical color 5 for cell labels or groups.",
    "Generic categorical color 6 for cell labels or groups."
  ),
  stringsAsFactors = FALSE
)
write_style_doc <- function(style_path) {
  doc_path <- sub("\\.tsv$", ".md", style_path)
  if (identical(doc_path, style_path)) doc_path <- paste0(style_path, ".md")
  lines <- c(
    "# SCENIC+ R Plot Style Parameters",
    "",
    "`plot_style_parameters.tsv` is the editable file. It has two columns: `parameter` and `value`.",
    "",
    "After editing the TSV, rerun all postprocess PDFs from the project root:",
    "",
    "```bash",
    "source project_env.sh",
    "spgrn-run-scenicplus-postprocess --task all --layer all",
    "```",
    "",
    "Supported rerun scope:",
    "",
    "- `--task audit`: output-tier audit only.",
    "- `--task figures`: regenerate source tables and 01-08 eRegulon figure PDFs.",
    "- `--task stats`: regenerate 09+ condition-statistics tables and PDFs.",
    "- `--task all`: complete postprocess output generation.",
    "- `--layer direct`, `--layer extended` or `--layer all`: choose the eRegulon layer; audit ignores this option.",
    "",
    "| parameter | default | meaning |",
    "|---|---:|---|"
  )
  rows <- merge(style_info, style_defaults, by = "parameter", sort = FALSE)
  rows <- rows[match(style_defaults$parameter, rows$parameter), ]
  rows$value <- gsub("\\|", "\\\\|", rows$value)
  rows$meaning <- gsub("\\|", "\\\\|", rows$meaning)
  lines <- c(lines, paste0("| `", rows$parameter, "` | `", rows$value, "` | ", rows$meaning, " |"))
  writeLines(lines, doc_path)
}
if (file.exists(style_path) && file.size(style_path) > 0) {
  sty_dt <- fread(style_path, sep = "\t", data.table = FALSE)
  sty <- setNames(style_defaults$value, style_defaults$parameter)
  if (all(c("parameter", "value") %in% names(sty_dt))) {
    valid <- sty_dt[as.character(sty_dt$parameter) %in% names(sty), , drop = FALSE]
    sty[as.character(valid$parameter)] <- as.character(valid$value)
  }
} else {
  sty <- setNames(style_defaults$value, style_defaults$parameter)
}
fwrite(data.table(parameter = names(sty), value = unname(sty)), style_path, sep = "\t")
write_style_doc(style_path)
sf <- function(k, default) suppressWarnings(as.numeric(ifelse(!is.na(sty[[k]]), sty[[k]], default)))
ss <- function(k, default) ifelse(!is.na(sty[[k]]), sty[[k]], default)
cap <- function(x, width = 105) paste(strwrap(x, width = width), collapse = "\n")

theme_set(
  theme_classic(base_size = sf("base_size", 6.5), base_family = ss("font_family", "Arial")) +
    theme(
      axis.text = element_text(size = sf("axis_text_size", 5.6), colour = "black"),
      plot.title = element_text(size = sf("plot_title_size", 7), face = "bold"),
      plot.caption = element_text(size = sf("caption_size", 5.0), hjust = 0),
      panel.grid = element_blank(),
      legend.frame = element_blank()
    )
)

save_pdf <- function(plot, path, width_mm = 120, height_mm = 90) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4, family = ss("font_family", "Arial"))
  } else {
    grDevices::pdf(path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Helvetica", useDingbats = FALSE)
  }
  print(plot)
  grDevices::dev.off()
  if (!file.exists(path) || file.size(path) == 0) stop("PDF was not created: ", path)
}

read_table <- function(stem, required = TRUE) {
  path <- file.path(outdir, named(stem, suffix, ".tsv"))
  if (!file.exists(path) || file.size(path) == 0) {
    if (required) stop("Missing table: ", path)
    return(NULL)
  }
  fread(path, sep = "\t", data.table = FALSE)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
condition_cols <- function(x, prefix) grep(paste0("^", prefix, "__"), names(x), value = TRUE)
condition_name <- function(x, prefix) sub(paste0("^", prefix, "__"), "", x)
read_priority_eregulons <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  path <- project_path(x)
  if (file.exists(path) && file.size(path) > 0) {
    tab <- fread(path, sep = "\t", data.table = FALSE)
    if (!any(c("eregulon", "display_label") %in% names(tab))) {
      stop("priority_eregulons file must contain `eregulon` or `display_label`: ", path)
    }
    return(tab)
  }
  vals <- trimws(unlist(strsplit(x, "[,;]")))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) return(NULL)
  data.frame(display_label = vals, stringsAsFactors = FALSE)
}
priority_labels <- function(tab, priority, label = NULL, layer = suffix) {
  if (is.null(priority) || nrow(tab) == 0) return(tab[0, , drop = FALSE])
  keep <- rep(TRUE, nrow(priority))
  if ("cell_label" %in% names(priority) && !is.null(label)) keep <- keep & priority$cell_label %in% label
  if ("layer" %in% names(priority) && nzchar(layer)) keep <- keep & priority$layer %in% layer
  pri <- priority[keep, , drop = FALSE]
  if (nrow(pri) == 0) return(tab[0, , drop = FALSE])
  hit <- rep(FALSE, nrow(tab))
  if ("eregulon" %in% names(pri) && "eregulon" %in% names(tab)) hit <- hit | tab$eregulon %in% pri$eregulon
  if ("display_label" %in% names(pri) && "display_label" %in% names(tab)) hit <- hit | tab$display_label %in% pri$display_label
  tab[hit, , drop = FALSE]
}
palette_conditions <- function(conditions) {
  pal <- setNames(rep(unname(sty[paste0("color_group_", 1:6)]), length.out = length(conditions)), conditions)
  ctrl_like <- grep("(^|[_-])(ctrl|control|wt|wildtype|wild_type)($|[_-])", conditions, ignore.case = TRUE, value = TRUE)
  if (length(ctrl_like) > 0) {
    pal[ctrl_like[1]] <- ss("color_control", "#2F4858")
    pal[setdiff(conditions, ctrl_like[1])[1]] <- ss("color_comparison", "#B57457")
  }
  pal
}

cluster_rows <- function(mat) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  if (nrow(mat) < 3) return(rownames(mat))
  row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
  if (all(row_sd == 0)) return(rownames(mat))
  rownames(mat)[stats::hclust(stats::dist(mat), method = "average")$order]
}

pdf_index <- 9
numbered <- function(stem) {
  path <- file.path(outdir, named(sprintf("%02d_%s", pdf_index, stem), suffix))
  pdf_index <<- pdf_index + 1
  path
}

sample_auc <- read_table("sample_mean_auc")
res <- read_table("condition_eregulon_auc_statistics")
priority <- read_priority_eregulons(args$priority_eregulons)

cond_col <- names(sample_auc)[2]
counts <- as.data.frame(table(sample_auc[[cond_col]]))
names(counts) <- c("condition", "n_samples")
p <- ggplot(counts, aes(condition, n_samples, fill = condition)) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = palette_conditions(counts$condition), guide = "none") +
  labs(x = NULL, y = "samples", title = "Samples per condition",
       caption = cap("Each bar is the number of biological samples contributing sample-level mean eRegulon AUC statistics."))
save_pdf(p, numbered("condition_sample_counts"), 88, 70)

if ("p_value" %in% names(res) && any(is.finite(res$p_value))) {
  res$neglog10_p <- -log10(pmax(as.numeric(res$p_value), .Machine$double.xmin))
  res$significant <- as.numeric(res$fdr) < 0.05
  plot_res <- res[is.finite(as.numeric(res$delta_mean_auc)) & is.finite(res$neglog10_p), ]
  lab <- unique(rbind(
    head(plot_res[order(plot_res$fdr, plot_res$p_value), ], 6),
    priority_labels(plot_res, priority)
  ))
  p <- ggplot(plot_res, aes(delta_mean_auc, neglog10_p)) +
    geom_point(aes(colour = significant), size = 0.8, alpha = 0.72) +
    geom_text(data = lab, aes(label = display_label), size = 1.8, vjust = -0.35, check_overlap = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "#B57457", `FALSE` = "#8F8F8F"), guide = "none") +
    labs(x = paste0("delta mean AUC (", unique(res$contrast)[1], ")"), y = "-log10 P", title = "Condition-level eRegulon AUC",
         caption = cap("Each point is one eRegulon tested on sample-level mean AUCell scores. Positive values follow the contrast shown on the x-axis."))
} else {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = "Insufficient sample replicates for P values") + theme_void()
}
save_pdf(p, numbered("condition_overall_eregulon_auc_volcano"), 120, 88)

top <- head(res$eregulon, min(25, nrow(res)))
mat_cols <- intersect(top, names(sample_auc))
if (length(mat_cols) > 0) {
  long <- data.table::melt(
    as.data.table(sample_auc[, c(names(sample_auc)[1:2], mat_cols), drop = FALSE]),
    id.vars = names(sample_auc)[1:2],
    variable.name = "eregulon",
    value.name = "AUC"
  )
  labels <- unique(res[, c("eregulon", "display_label")])
  long <- merge(long, labels, by = "eregulon", all.x = TRUE)
  long$.sample_for_cluster <- long[[names(sample_auc)[1]]]
  wide <- xtabs(AUC ~ eregulon + .sample_for_cluster, data = long)
  row_order <- cluster_rows(wide)
  long$sample_condition <- paste(long[[names(sample_auc)[1]]], long[[cond_col]], sep = "|")
  long$display_label <- factor(long$display_label, levels = rev(unique(labels$display_label[match(row_order, labels$eregulon)])))
  p <- ggplot(long, aes(sample_condition, display_label, fill = AUC)) +
    geom_tile() +
    scale_fill_gradient(low = "#F2F4F5", high = "#2F4858") +
    labs(x = NULL, y = NULL, title = "Top eRegulons: sample mean AUC",
         caption = cap("Rows are top condition-associated eRegulons and columns are biological samples annotated by condition.")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_pdf(p, numbered("condition_overall_top_sample_mean_auc_heatmap"), 170, 128)
}

by_label <- read_table("by_cell_label_condition_eregulon_auc_statistics", required = FALSE)
if (!is.null(by_label) && nrow(by_label) > 0) {
  label_col <- names(by_label)[1]
  by_label$abs_delta <- abs(as.numeric(by_label$delta_mean_auc))
  top_regs <- unique(by_label$eregulon[order(-by_label$abs_delta)])[1:min(30, length(unique(by_label$eregulon)))]
  hm <- by_label[by_label$eregulon %in% top_regs, ]
  labels <- unique(hm[, c("eregulon", "display_label")])
  hm$.label_for_cluster <- hm[[label_col]]
  wide <- xtabs(as.numeric(delta_mean_auc) ~ eregulon + .label_for_cluster, data = hm)
  row_order <- cluster_rows(wide)
  hm$display_label <- factor(hm$display_label, levels = rev(unique(labels$display_label[match(row_order, labels$eregulon)])))
  p <- ggplot(hm, aes(.data[[label_col]], display_label, fill = as.numeric(delta_mean_auc))) +
    geom_tile() +
    scale_fill_gradient2(low = "#355C8C", mid = "white", high = "#B2574E", midpoint = 0, name = "delta AUC") +
    labs(x = NULL, y = NULL, title = "Condition effect by cell label",
         caption = cap("Each tile is the condition effect size for one eRegulon within one cell label, computed from sample-level mean AUCell scores.")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p, numbered("condition_by_cell_label_effect_heatmap"), 138, 128)
  for (label in unique(by_label[[label_col]])) {
    sub <- by_label[by_label[[label_col]] == label, ]
    if ("p_value" %in% names(sub) && any(is.finite(sub$p_value))) {
      sub$neglog10_p <- -log10(pmax(as.numeric(sub$p_value), .Machine$double.xmin))
      sub$significant <- as.numeric(sub$fdr) < 0.05
      plot_sub <- sub[is.finite(as.numeric(sub$delta_mean_auc)) & is.finite(sub$neglog10_p), ]
      lab <- unique(rbind(
        head(plot_sub[order(plot_sub$fdr, plot_sub$p_value), ], 5),
        priority_labels(plot_sub, priority, label = label)
      ))
      p <- ggplot(plot_sub, aes(delta_mean_auc, neglog10_p)) +
        geom_point(aes(colour = significant), size = 0.8, alpha = 0.72) +
        geom_text(data = lab, aes(label = display_label), size = 1.8, vjust = -0.35, check_overlap = FALSE) +
        scale_colour_manual(values = c(`TRUE` = "#B57457", `FALSE` = "#8F8F8F"), guide = "none") +
        labs(x = paste0("delta mean AUC (", unique(sub$contrast)[1], ")"), y = "-log10 P", title = paste0(label, ": condition eRegulon AUC"),
             caption = cap("Each point is one eRegulon tested within the displayed cell label using sample-level mean AUCell scores."))
      save_pdf(p, numbered(paste0("condition_by_cell_label_eregulon_auc_volcano__", safe_name(label))), 120, 88)
    }
  }
}

cat("SCENIC+ R-rendered condition statistics figures written to", outdir, "\n")
