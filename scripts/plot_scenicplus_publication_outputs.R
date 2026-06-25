#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  out <- list(
    outdir = "results/scenicplus_figures",
    file_suffix = "",
    plot_style_config = "results/scenicplus_figures/plot_style_parameters.tsv"
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
project_path <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return("")
  if (grepl("^/", x)) normalizePath(x, mustWork = FALSE) else normalizePath(file.path(project_dir, x), mustWork = FALSE)
}

named <- function(stem, suffix, ext = ".pdf") {
  if (is.null(suffix) || suffix == "") paste0(stem, ext) else paste0(stem, "_", suffix, ext)
}

style_defaults <- data.table(
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
  )
)

style_info <- data.table(
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
  )
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

load_style <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  defaults <- setNames(style_defaults$value, style_defaults$parameter)
  if (!file.exists(path) || file.size(path) == 0) {
    fwrite(style_defaults, path, sep = "\t")
    write_style_doc(path)
    return(defaults)
  }
  x <- fread(path, sep = "\t", data.table = FALSE)
  if (!all(c("parameter", "value") %in% names(x))) {
    fwrite(style_defaults, path, sep = "\t")
    write_style_doc(path)
    return(defaults)
  }
  valid <- x[as.character(x$parameter) %in% names(defaults), , drop = FALSE]
  defaults[as.character(valid$parameter)] <- as.character(valid$value)
  legacy <- c(
    color_group_1 = "color_fallback_1",
    color_group_2 = "color_fallback_2",
    color_group_3 = "color_fallback_3",
    color_group_4 = "color_fallback_4",
    color_group_5 = "color_fallback_5",
    color_group_6 = "color_fallback_6"
  )
  old <- setNames(as.character(x$value), as.character(x$parameter))
  for (new_name in names(legacy)) {
    if (legacy[[new_name]] %in% names(old) && !(new_name %in% as.character(x$parameter))) {
      defaults[[new_name]] <- old[[legacy[[new_name]]]]
    }
  }
  fwrite(data.table(parameter = names(defaults), value = unname(defaults)), path, sep = "\t")
  write_style_doc(path)
  defaults
}

args <- parse_args()
outdir <- project_path(args$outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
style_path <- project_path(args$plot_style_config)
sty <- load_style(style_path)

sf <- function(k, default) suppressWarnings(as.numeric(ifelse(!is.na(sty[[k]]), sty[[k]], default)))
si <- function(k, default) suppressWarnings(as.integer(round(as.numeric(ifelse(!is.na(sty[[k]]), sty[[k]], default)))))
ss <- function(k, default) ifelse(!is.na(sty[[k]]), sty[[k]], default)
suffix <- args$file_suffix
cap <- function(x, width = 115) paste(strwrap(x, width = width), collapse = "\n")

theme_set(
  theme_classic(base_size = sf("base_size", 6.5), base_family = ss("font_family", "Arial")) +
    theme(
      axis.line = element_line(linewidth = sf("axis_line_width", 0.30), colour = "black"),
      axis.ticks = element_line(linewidth = sf("axis_line_width", 0.30), colour = "black"),
      axis.text = element_text(size = sf("axis_text_size", 5.6), colour = "black"),
      axis.title = element_text(size = sf("axis_title_size", 6.3), colour = "black"),
      legend.title = element_text(size = sf("legend_title_size", 5.8)),
      legend.text = element_text(size = sf("legend_text_size", 5.5)),
      strip.text = element_text(size = sf("strip_text_size", 6.2), face = "bold"),
      plot.title = element_text(size = sf("plot_title_size", 7), face = "bold"),
      plot.caption = element_text(size = sf("caption_size", 5.0), hjust = 0, colour = "black"),
      panel.grid = element_blank()
    )
)

save_pdf <- function(plot, path, width_mm = 183, height_mm = 120) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width, height = height, family = ss("font_family", "Arial"))
  } else {
    grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  }
  print(plot)
  grDevices::dev.off()
  if (!file.exists(path) || file.size(path) == 0) stop("PDF was not created: ", path)
}

read_source <- function(stem, required = TRUE) {
  path <- file.path(outdir, named(stem, suffix, ".tsv"))
  if (!file.exists(path) || file.size(path) == 0) {
    if (required) stop("Missing source table: ", path)
    return(NULL)
  }
  fread(path, sep = "\t", data.table = FALSE)
}

regulon_levels <- function(x) unique(as.character(x$display_label[order(x$eregulon)]))
group_palette <- function(groups) {
  colors <- unname(sty[paste0("color_group_", 1:6)])
  setNames(rep(colors, length.out = length(groups)), groups)
}

add_plot_labels <- function(src) {
  labs <- unique(src[, c("eregulon", "display_label")])
  labs$plot_label <- make.unique(as.character(labs$display_label))
  merge(src, labs[, c("eregulon", "plot_label")], by = "eregulon", all.x = TRUE, sort = FALSE)
}

cluster_levels <- function(src, row_col, col_col, value_col, cluster_cols = TRUE) {
  mat <- xtabs(as.numeric(src[[value_col]]) ~ src[[row_col]] + src[[col_col]])
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  row_levels <- rownames(mat)
  col_levels <- colnames(mat)
  if (nrow(mat) >= 3 && any(apply(mat, 1, stats::sd) > 0)) {
    row_levels <- rownames(mat)[stats::hclust(stats::dist(mat), method = "average")$order]
  }
  if (cluster_cols && ncol(mat) >= 3 && any(apply(mat, 2, stats::sd) > 0)) {
    col_levels <- colnames(mat)[stats::hclust(stats::dist(t(mat)), method = "average")$order]
  }
  list(rows = row_levels, cols = col_levels)
}

plot_heatmap <- function(src, value_col, title, legend, path, diverging = TRUE, fixed_col_order = FALSE) {
  src <- add_plot_labels(src)
  ord <- cluster_levels(src, "plot_label", "group", value_col, cluster_cols = !fixed_col_order)
  src$plot_label <- factor(src$plot_label, levels = rev(ord$rows))
  if (fixed_col_order) {
    src$group <- factor(src$group, levels = unique(src$group))
  } else {
    src$group <- factor(src$group, levels = ord$cols)
  }
  p <- ggplot(src, aes(x = group, y = plot_label, fill = .data[[value_col]])) +
    geom_tile(linewidth = sf("heatmap_tile_border_width", 0), colour = NA) +
    labs(x = NULL, y = NULL, title = title, caption = cap(legend), fill = value_col) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = sf("axis_text_size", 5.6))
    )
  if (diverging) {
    vmax <- max(abs(src[[value_col]]), na.rm = TRUE)
    if (!is.finite(vmax) || vmax == 0) vmax <- 1
    p <- p + scale_fill_gradient2(low = "#355C8C", mid = "white", high = "#B2574E", midpoint = 0, limits = c(-vmax, vmax), oob = scales::squish)
  } else {
    p <- p + scale_fill_gradient(low = "#F2F4F5", high = "#2F4858")
  }
  save_pdf(p, path, width_mm = ifelse(fixed_col_order, 195, 160), height_mm = 135)
}

plot_dot <- function(src, title, legend, path) {
  fixed_col_order <- grepl("_condition_", basename(path))
  src <- add_plot_labels(src)
  ord <- cluster_levels(src, "plot_label", "group", "z_auc", cluster_cols = !fixed_col_order)
  src$plot_label <- factor(src$plot_label, levels = ord$rows)
  if (fixed_col_order) {
    src$group <- factor(src$group, levels = unique(src$group[order(src$group_order)]))
  } else {
    src$group <- factor(src$group, levels = rev(ord$cols))
  }
  vmax <- max(abs(src$z_auc), na.rm = TRUE)
  if (!is.finite(vmax) || vmax == 0) vmax <- 1
  p <- ggplot(src, aes(x = plot_label, y = group)) +
    geom_point(aes(size = active_fraction, colour = z_auc), alpha = 0.92) +
    scale_size(range = c(sf("dot_min_size", 0.7), sf("dot_size_scale", 3.0)), name = "active fraction") +
    scale_colour_gradient2(low = "#355C8C", mid = "white", high = "#B2574E", midpoint = 0, limits = c(-vmax, vmax), oob = scales::squish, name = "AUC z-score") +
    labs(x = NULL, y = NULL, title = title, caption = cap(legend)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_pdf(p, path, width_mm = 195, height_mm = 118)
}

plot_auc_umap <- function(src, path) {
  labels <- unique(src$display_label)
  chunks <- split(labels, ceiling(seq_along(labels) / 6))
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = 183 / 25.4, height = 118 / 25.4, family = ss("font_family", "Arial"))
  } else {
    grDevices::pdf(path, width = 183 / 25.4, height = 118 / 25.4, family = "Helvetica", useDingbats = FALSE)
  }
  for (chunk in chunks) {
    x <- src[src$display_label %in% chunk, , drop = FALSE]
    p <- ggplot(x, aes(umap_1, umap_2, colour = auc)) +
      geom_point(size = sf("umap_point_size", 0.20), alpha = sf("umap_point_alpha", 0.75), stroke = 0) +
      facet_wrap(~display_label, ncol = 3) +
      scale_colour_gradient(low = "#F0E6E1", high = "#B57457") +
      labs(x = NULL, y = NULL, title = "eRegulon AUCell activity on annotated embedding",
           caption = cap("Each panel projects single-cell eRegulon activity onto the annotated embedding. Color is AUCell activity for the displayed eRegulon.")) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    print(p)
  }
  grDevices::dev.off()
}

plot_auc_umap_condition <- function(src, path) {
  labels <- unique(src$display_label)
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = 183 / 25.4, height = 118 / 25.4, family = ss("font_family", "Arial"))
  } else {
    grDevices::pdf(path, width = 183 / 25.4, height = 118 / 25.4, family = "Helvetica", useDingbats = FALSE)
  }
  for (lab in labels) {
    x <- src[src$display_label == lab, , drop = FALSE]
    p <- ggplot(x, aes(umap_1, umap_2, colour = auc)) +
      geom_point(size = sf("umap_point_size", 0.20), alpha = sf("umap_point_alpha", 0.75), stroke = 0) +
      facet_wrap(~condition, ncol = min(si("max_panels_per_row", 3), length(unique(x$condition)))) +
      scale_colour_gradient(low = "#F0E6E1", high = "#B57457") +
      labs(x = NULL, y = NULL, title = paste0(lab, ": condition-resolved AUCell activity"),
           caption = cap("Condition panels use the same embedding and color scale for the displayed eRegulon; each panel contains cells from one condition.")) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    print(p)
  }
  grDevices::dev.off()
}

plot_activity_embedding <- function(src, path) {
  groups <- unique(src$cell_label)
  pal <- group_palette(groups)
  conds <- unique(src$condition[src$condition != ""])
  make_panel <- function(label) {
    if (label == "All") {
      ggplot(src, aes(eregulon_umap_1, eregulon_umap_2, colour = cell_label)) +
        geom_point(size = sf("embedding_point_size", 0.25), alpha = sf("embedding_point_alpha", 0.75), stroke = 0) +
        scale_colour_manual(values = pal, drop = FALSE) +
        labs(title = "All cells", x = "eRegulon UMAP 1", y = "eRegulon UMAP 2", colour = NULL)
    } else {
      fg <- src[src$condition == label, , drop = FALSE]
      ggplot() +
        geom_point(data = src, aes(eregulon_umap_1, eregulon_umap_2), colour = "#D0D0D0", size = sf("condition_background_point_size", 0.10), alpha = sf("condition_background_alpha", 0.18), stroke = 0) +
        geom_point(data = fg, aes(eregulon_umap_1, eregulon_umap_2, colour = cell_label), size = sf("embedding_point_size", 0.25), alpha = sf("embedding_point_alpha", 0.75), stroke = 0) +
        scale_colour_manual(values = pal, drop = FALSE) +
        labs(title = label, x = "eRegulon UMAP 1", y = "eRegulon UMAP 2", colour = NULL)
    }
  }
  panels <- c("All", conds)
  plots <- lapply(panels, make_panel)
  p <- wrap_plots(plots, ncol = min(si("max_panels_per_row", 3), length(plots)), guides = "collect") +
    plot_annotation(
      title = "eRegulon activity embedding",
      caption = cap("The embedding is computed from the eRegulon AUCell matrix. The first panel shows all cells; condition panels reuse the same coordinates and display each condition separately.")
    )
  save_pdf(p, path, width_mm = 183, height_mm = 68 * ceiling(length(plots) / min(si("max_panels_per_row", 3), length(plots))))
}

plot_region_gene <- function(counts, dist_imp, path) {
  p1 <- ggplot(counts, aes(n_regions)) +
    geom_histogram(binwidth = 1, fill = "#667A9A", colour = "white", linewidth = 0.1) +
    labs(x = "regions per gene", y = "genes", title = "Enhancer load per gene")
  if (!is.null(dist_imp) && nrow(dist_imp) > 0) {
    dist_imp$log_distance <- log10(dist_imp$abs_distance + 1)
    p2 <- ggplot(dist_imp, aes(log_distance, importance)) +
      geom_bin2d(bins = 42) +
      scale_fill_gradient(low = "#F2F4F5", high = "#2F4858") +
      labs(x = "log10 distance to gene + 1", y = "region-gene importance", title = "Distance vs link score")
  } else {
    p2 <- ggplot() + annotate("text", x = 0, y = 0, label = "distance or importance column not found") + theme_void()
  }
  p <- p1 + p2 + plot_annotation(caption = cap("Region-to-gene links are model-level SCENIC+ associations. The panels summarize enhancer load and distance/importance structure; condition effects are evaluated with AUCell statistics."))
  save_pdf(p, path, width_mm = 183, height_mm = 80)
}

plot_overlap <- function(src, path) {
  wide <- xtabs(jaccard ~ display_a + display_b, data = src)
  mat <- as.matrix(wide)
  mat[!is.finite(mat)] <- 0
  if (nrow(mat) >= 3 && any(apply(mat, 1, stats::sd) > 0)) {
    ord <- rownames(mat)[stats::hclust(stats::dist(mat), method = "average")$order]
  } else {
    ord <- rownames(mat)
  }
  src$display_a <- factor(src$display_a, levels = rev(ord))
  src$display_b <- factor(src$display_b, levels = ord)
  p <- ggplot(src, aes(display_b, display_a, fill = jaccard)) +
    geom_tile() +
    scale_fill_gradient(low = "#F4ECEA", high = "#8E5A63", name = "Jaccard") +
    labs(x = NULL, y = NULL, title = "Target-region overlap",
         caption = cap("Each cell is the Jaccard overlap between two eRegulons' linked regulatory-region sets. This is a model-level structural view.")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_pdf(p, path, width_mm = 160, height_mm = 150)
}

plot_network <- function(edges, path) {
  if (is.null(edges) || nrow(edges) == 0) {
    p <- ggplot() + annotate("text", x = 0, y = 0, label = "No TF-target edges found") + theme_void()
    save_pdf(p, path, width_mm = 120, height_mm = 80)
    return()
  }
  tfs <- unique(edges$tf)
  tg <- unique(edges$target_gene)
  nodes <- rbind(
    data.frame(node = tfs, type = "TF", x = 0, y = seq_along(tfs)),
    data.frame(node = tg, type = "target", x = 1, y = seq_along(tg) / max(1, length(tg)) * max(1, length(tfs)))
  )
  e <- merge(edges, nodes[, c("node", "x", "y")], by.x = "tf", by.y = "node")
  names(e)[names(e) %in% c("x", "y")] <- c("x_tf", "y_tf")
  e <- merge(e, nodes[, c("node", "x", "y")], by.x = "target_gene", by.y = "node")
  names(e)[names(e) %in% c("x", "y")] <- c("x_target", "y_target")
  p <- ggplot() +
    geom_segment(data = e, aes(x = x_tf, y = y_tf, xend = x_target, yend = y_target, linewidth = weight), colour = "#777777", alpha = 0.25) +
    geom_point(data = nodes, aes(x, y, fill = type, size = type), shape = 21, colour = "white", stroke = 0.15) +
    geom_text(data = nodes[nodes$type == "TF", ], aes(x - 0.02, y, label = node), hjust = 1, size = 2) +
    geom_text(data = nodes[nodes$type == "target", ], aes(x + 0.02, y, label = node), hjust = 0, size = 1.35, check_overlap = TRUE) +
    scale_fill_manual(values = c(TF = "#B57457", target = "#9AA6B2")) +
    scale_size_manual(values = c(TF = 2.7, target = 1.1)) +
    scale_linewidth(range = c(0.15, 0.8), guide = "none") +
    coord_cartesian(xlim = c(-0.18, 1.42), clip = "off") +
    labs(title = "TF-target network", caption = cap("The network displays selected TF-target links. TF nodes are shown at left; target genes are shown at right; edge width follows the available link weight.")) +
    theme_void(base_size = sf("base_size", 6.5)) +
    theme(legend.position = "none", plot.caption = element_text(size = sf("caption_size", 5.0), hjust = 0), plot.margin = margin(5.5, 45, 5.5, 18))
  save_pdf(p, path, width_mm = 210, height_mm = 135)
}

auc <- read_source("source_auc_heatmap_zscore")
plot_heatmap(auc, "z_auc", "eRegulon activity", "Rows are eRegulons and columns are cell labels. Values are row-scaled AUCell means.", file.path(outdir, named("01_eregulon_auc_heatmap", suffix)), TRUE, FALSE)

rss <- read_source("source_rss_specificity")
plot_heatmap(rss, "specificity", "eRegulon specificity", "Rows are eRegulons and columns are cell labels. Values summarize regulon specificity by cell label.", file.path(outdir, named("02_eregulon_specificity_heatmap", suffix)), FALSE, FALSE)

dot <- read_source("source_dot_heatmap")
plot_dot(dot, "eRegulon activity dot heatmap", "Dot color shows row-scaled mean AUCell activity and dot size shows the fraction of active cells.", file.path(outdir, named("03_eregulon_dot_heatmap", suffix)))

auc_cond <- read_source("source_auc_heatmap_zscore_condition", required = FALSE)
if (!is.null(auc_cond)) {
  plot_heatmap(auc_cond, "z_auc", "eRegulon activity by state and condition", "Columns are ordered as cell label followed by condition, allowing condition effects to be compared within each cell state.", file.path(outdir, named("01_eregulon_auc_heatmap", paste0("condition_", suffix))), TRUE, TRUE)
}
rss_cond <- read_source("source_rss_specificity_condition", required = FALSE)
if (!is.null(rss_cond)) {
  plot_heatmap(rss_cond, "specificity", "eRegulon specificity by state and condition", "Specificity is calculated against cell-state/condition strata; columns preserve cell label then condition ordering.", file.path(outdir, named("02_eregulon_specificity_heatmap", paste0("condition_", suffix))), FALSE, TRUE)
}
dot_cond <- read_source("source_dot_heatmap_condition", required = FALSE)
if (!is.null(dot_cond)) {
  plot_dot(dot_cond, "eRegulon activity by state and condition", "Each dot summarizes one eRegulon in one cell-state/condition stratum.", file.path(outdir, named("03_eregulon_dot_heatmap", paste0("condition_", suffix))))
}

umap <- read_source("source_eregulon_auc_umap", required = FALSE)
if (!is.null(umap)) {
  plot_auc_umap(umap, file.path(outdir, named("04_eregulon_auc_umap", suffix)))
  if ("condition" %in% names(umap) && length(unique(umap$condition[umap$condition != ""])) > 1) {
    plot_auc_umap_condition(umap, file.path(outdir, named("04_eregulon_auc_umap", paste0("condition_", suffix))))
  }
}

emb <- read_source("source_eregulon_activity_embedding")
plot_activity_embedding(emb, file.path(outdir, named("05_eregulon_activity_embedding", suffix)))

counts <- read_source("source_regions_per_gene", required = FALSE)
dist_imp <- read_source("source_region_distance_importance", required = FALSE)
if (!is.null(counts)) plot_region_gene(counts, dist_imp, file.path(outdir, named("06_region_gene_link_structure", suffix)))

overlap <- read_source("source_target_region_overlap", required = FALSE)
if (!is.null(overlap)) plot_overlap(overlap, file.path(outdir, named("07_target_region_overlap_heatmap", suffix)))

edges <- read_source("source_tf_target_network_edges", required = FALSE)
plot_network(edges, file.path(outdir, named("08_tf_target_network", suffix)))

cat("SCENIC+ R-rendered figures written to", outdir, "\n")
