
# =======================================
# Volcano Plot (EdgeR)
# =======================================
# Construct a volcano plot from a
# differential expression results table.
plot_volcano_edge_r <- function(
    de_table,
    title = NULL,
    fdr_cutoff = 0.05,
    logfc_cutoff = 1,
    n_labels = 10
) {
  
  df <- de_table

  # Classify genes
  
  df$Significance <- "Not significant"
  
  df$Significance[
    df$FDR < fdr_cutoff &
      df$logFC > 0
  ] <- "Up"
  
  df$Significance[
    df$FDR < fdr_cutoff &
      df$logFC < 0
  ] <- "Down"
  
  # Top genes to label
  
  sig <- subset(
    df,
    FDR < fdr_cutoff &
      abs(logFC) > logfc_cutoff
  )
  
  top_sig <- NULL
  
  if (nrow(sig) > 0) {
    
    top_sig <- sig[
      order(sig$FDR),
      ,
      drop = FALSE
    ]
    
    top_sig <- head(
      top_sig,
      n_labels
    )
    
    top_sig$Gene <- rownames(top_sig)
    
  }
  
  # Volcano plot
  
  p <-
    
    ggplot(
      df,
      aes(
        x = logFC,
        y = -log10(FDR),
        colour = Significance
      )
    ) +
    
    geom_point(
      alpha = 0.7,
      size = 1.5
    ) +
    
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      colour = "grey60"
    ) +
    
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      colour = "grey60"
    ) +
    
    scale_colour_manual(
      values = c(
        Up = "#D55E00",
        Down = "#0072B2",
        `Not significant` = "grey80"
      )
    ) +
    
    labs(
      title = title,
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P)
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  # Add gene labels
  if (!is.null(top_sig)) {
    
    p <- p +
      
      ggrepel::geom_text_repel(
        data = top_sig,
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "grey60",
        show.legend = FALSE
      )
    
  }
  
  p
  
}
  
# =======================================
# GO Dotplot (EdgeR)
# =======================================
# Construct a Gene Ontology enrichment
# dotplot.

plot_go_edgeR <- function(
    go_result,
    title = NULL,
    n_terms = 10
) {
  
  if (is.null(go_result))
    return(NULL)
  
  p <-
    enrichplot::dotplot(
      go_result,
      showCategory = n_terms
    ) +
    
    scale_y_discrete(
      labels = \(x)
      stringr::str_wrap(
        x,
        width = 30
      )
    ) +
    
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.text.x = element_text(
        colour = "black"
      ),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "black",
        linewidth = 0.6
      ),
      legend.title = element_text(
        face = "bold"
      ),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
  
  return(p)
}

# =======================================
# Volcano Plot (Seurat)
# =======================================

plot_volcano_seurat <- function(
    de_table,
    title = NULL,
    fdr_cutoff = 0.05,
    logfc_cutoff = 1,
    n_labels = 10,
    ymax = 50
) {
  
  df <- de_table
  
  df$Significance <- "Not significant"
  
  df$Significance[
    df$p_val_adj < fdr_cutoff &
      df$avg_log2FC > logfc_cutoff
  ] <- "Up"
  
  df$Significance[
    df$p_val_adj < fdr_cutoff &
      df$avg_log2FC < -logfc_cutoff
  ] <- "Down"
  
  df$neglog10_raw <- -log10(df$p_val_adj)
  
  df$clipped <- is.infinite(df$neglog10_raw) | df$neglog10_raw > ymax

  df$neglog10 <- pmin(df$neglog10_raw, ymax)
  
  n_clipped <- sum(df$clipped, na.rm = TRUE)
  if (n_clipped > 0) {
    message(sprintf(
      "plot_volcano_seurat: %d genes have p_val_adj underflowed to 0 
        or -log10(p) > ymax (%d): capped at y = %d and shown as triangles.",
      n_clipped, ymax, ymax
    ))
  }
  
  sig <- subset(
    df,
    p_val_adj < fdr_cutoff &
      abs(avg_log2FC) > logfc_cutoff
  )
  
  top_sig <- NULL
  
  if (nrow(sig) > 0) {
    
    top_sig <- sig[order(-abs(sig$avg_log2FC), sig$p_val_adj), ]
    top_sig <- head(top_sig, n_labels)
    top_sig$Gene <- rownames(top_sig)
    top_sig$neglog10 <- pmin(top_sig$neglog10_raw, ymax)
    
  }
  
  p <- ggplot(
    df,
    aes(
      x = avg_log2FC,
      y = neglog10,
      colour = Significance,
      shape = clipped
    )
  ) +
    
    geom_point(alpha = 0.7, size = 1.5) +
    
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed"
    ) +
    
    coord_cartesian(
      xlim = c(-3, 3),
      ylim = c(0, ymax)
    ) +
    
    scale_colour_manual(
      values = c(
        Up = "#D55E00",
        Down = "#0072B2",
        `Not significant` = "grey80"
      )
    ) +
    
    scale_shape_manual(
      values = c('FALSE' = 16, 'TRUE' = 17),
      guide = "none"
    ) +
    
    labs(
      title = title,
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~(capped~at.(ymax)))
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  if (!is.null(top_sig)) {
    
    p <- p +
      ggrepel::geom_text_repel(
        data = top_sig,
        aes(label = Gene),
        max.overlaps = Inf
      )
    
  }
  
  p
}

# ===================================
# GO Dotplot (Seurat)
# ===================================
plot_go_seurat <- function(
    go_result,
    title = NULL,
    n_terms = 15
){
  
  if (is.null(go_result))
    return(NULL)
  
  enrichplot::dotplot(
    go_result,
    showCategory = n_terms
  ) +
    
    scale_y_discrete(
      labels = \(x)
      stringr::str_wrap(
        x,
        width = 25
      )
    ) +
    
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text.y = element_text(
        colour = "black"
      ),
      
      axis.text.x = element_text(
        colour = "black"
      ),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        linewidth = 0.6
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.background = element_blank(),
      
      legend.key = element_blank()
      
    )
  
}

# =======================================
# Top DEGs Bar Chart (Seurat)
# =======================================
# Ranks genes by effect size (avg_log2FC) rather than p-value, since
# per-cell Wilcoxon tests on large Xenium cell counts saturate p_val_adj.

plot_top_genes_bar <- function(
    de_table,
    title = NULL,
    run_name = NULL,
    fdr_cutoff = 0.05,
    logfc_cutoff = 1,
    n_genes = 20,
    min_pct = 0.1
) {
  
  if (!is.data.frame(de_table)) {
    stop(
      "plot_top_genes_bar: `de_table` must be a single data.frame ",
      "(e.g. de_results[[\"comparison_name\"]]), not a list of tables. ",
      "Got class: ", paste(class(de_table), collapse = ", ")
    )
  }
  
  if (!all(c("pct.1", "pct.2") %in% colnames(de_table))) {
    stop(
      "plot_top_genes_bar: `de_table` must contain pct.1 and pct.2 columns ",
      "(standard Seurat FindMarkers output). Found columns: ",
      paste(colnames(de_table), collapse = ", ")
    )
  }
  
  df <- de_table
  df$Gene <- rownames(df)
  
  df$Significance <- "Not significant"
  df$Significance[
    df$p_val_adj < fdr_cutoff & df$avg_log2FC > logfc_cutoff
  ] <- "Up"
  df$Significance[
    df$p_val_adj < fdr_cutoff & df$avg_log2FC < -logfc_cutoff
  ] <- "Down"
  
  ## Flag genes whose fold change is likely inflated by a near-zero
  ## pseudocount denominator rather than reflecting a robust shift.
  df$low_pct <- (df$pct.1 < min_pct) & (df$pct.2 < min_pct)
  
  sig <- subset(df, Significance %in% c("Up", "Down"))
  
  n_low_pct_excluded <- sum(sig$low_pct)
  
  sig <- subset(sig, !low_pct)
  
  if (nrow(sig) == 0) {
    warning(
      "plot_top_genes_bar: no genes pass fdr_cutoff/logfc_cutoff/min_pct; ",
      "returning NULL."
    )
    return(NULL)
  }
  
  n_half <- ceiling(n_genes / 2)
  
  up <- sig[sig$Significance == "Up", ]
  up <- up[order(-up$avg_log2FC), ]
  up <- head(up, n_half)
  
  down <- sig[sig$Significance == "Down", ]
  down <- down[order(down$avg_log2FC), ]
  down <- head(down, n_half)
  
  top <- rbind(up, down)
  
  ## Order genes on the y-axis by effect size
  top$Gene <- factor(top$Gene, levels = top$Gene[order(top$avg_log2FC)])
  
  ## Label each bar with pct.1 / pct.2 to show underlying detection rates.
  top$pct_label <- sprintf("%.0f%% / %.0f%%", top$pct.1 * 100, top$pct.2 * 100)
  
  n_up_total <- sum(sig$Significance == "Up")
  n_down_total <- sum(sig$Significance == "Down")
  
  if (!is.null(run_name)) {
    title <- paste(
      gsub("_", " ", run_name),
      gsub("_", " ", title),
      sep = " — "
    )
  }
  
  p <- ggplot(
    top,
    aes(x = avg_log2FC, y = Gene, fill = Significance)
  ) +
    
    geom_col(width = 0.7) +
    
    geom_vline(xintercept = 0, colour = "black") +
    
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      colour = "grey40"
    ) +
    
    geom_text(
      aes(
        label = pct_label,
        hjust = ifelse(avg_log2FC > 0, -0.05, 1.05)
      ),
      size = 3,
      colour = "grey30"
    ) +
    
    scale_fill_manual(
      values = c(Up = "#D55E00", Down = "#0072B2"),
      guide = "none"
    ) +
    
    ## Give the pct labels room to sit past the longest bars without
    ## being clipped by the panel edge.
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    
    labs(
      title = title,
      subtitle = sprintf(
        paste0(
          "%d up / %d down total (p_val_adj < %s, |log2FC| > %s) -- top %d shown\n",
          "%d gene(s) excluded: pct.1 & pct.2 both < %s%% (likely pseudocount-inflated fold change).\n ",
          "Bar labels show pct.1 / pct.2 (%% cells expressing)."
        ),
        n_up_total, n_down_total, fdr_cutoff, logfc_cutoff, nrow(top),
        n_low_pct_excluded, min_pct * 100
      ),
      x = expression(log[2]~Fold~Change),
      y = NULL
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30"),
      axis.title.x = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(size = 9)
    )
  
  p
}

# ======================================================
# DE Gene Heatmap Across Time Points (Seurat)
# ======================================================
# Uses the same significance/effect-size/pct filtering logic as
# plot_top_genes_bar() to build the gene set, so the two plots stay
# consistent with each other.

plot_de_heatmap <- function(
    de_results,             
    avg_expr = NULL,             # OPTIONAL: precomputed genes x timepoints
    seurat_obj = NULL,          
    group_by = NULL,               
    assay = NULL,                   
    slot = "data",                   
    timepoint_order,             # e.g. c("1wk", "2wk", "4wk", "12wk")
    title = NULL,
    fdr_cutoff = 0.05,
    logfc_cutoff = 1,
    min_pct = 0.1,
    n_genes_per_comparison = 20  # top up+down per comparison, before union
) {
  

  if (!is.list(de_results) || is.data.frame(de_results)) {
    stop("plot_de_heatmap: `de_results` must be a named list of per-comparison data.frames.")
  }
  
  if (!all(vapply(de_results, is.data.frame, logical(1)))) {
    stop("plot_de_heatmap: every element of `de_results` must be a data.frame.")
  }
  
  if (!all(c("pct.1", "pct.2") %in% unlist(lapply(de_results, colnames)))) {
    stop("plot_de_heatmap: every de_results table must contain pct.1 and pct.2 columns.")
  }
  
  if (is.null(avg_expr) && (is.null(seurat_obj) || is.null(group_by))) {
    stop(
      "plot_de_heatmap: supply either `avg_expr` directly, or both ",
      "`seurat_obj` and `group_by` so it can be computed internally."
    )
  }
  
  get_top_genes <- function(df, comparison_name) {
    
    df$Gene <- rownames(df)
    df$Significance <- "Not significant"
    df$Significance[df$p_val_adj < fdr_cutoff & df$avg_log2FC > logfc_cutoff] <- "Up"
    df$Significance[df$p_val_adj < fdr_cutoff & df$avg_log2FC < -logfc_cutoff] <- "Down"
    df$low_pct <- (df$pct.1 < min_pct) & (df$pct.2 < min_pct)
    
    sig <- subset(df, Significance %in% c("Up", "Down") & !low_pct)
    
    if (nrow(sig) == 0) return(NULL)
    
    n_half <- ceiling(n_genes_per_comparison / 2)
    
    up <- head(sig[sig$Significance == "Up", ][order(-sig[sig$Significance == "Up", "avg_log2FC"]), ], n_half)
    down <- head(sig[sig$Significance == "Down", ][order(sig[sig$Significance == "Down", "avg_log2FC"]), ], n_half)
    
    genes <- c(up$Gene, down$Gene)
    if (length(genes) == 0) return(NULL)
    
    data.frame(Gene = genes, source = comparison_name, stringsAsFactors = FALSE)
  }
  
  gene_hits <- Map(get_top_genes, de_results, names(de_results))
  gene_hits <- do.call(rbind, gene_hits[!vapply(gene_hits, is.null, logical(1))])
  
  if (is.null(gene_hits) || nrow(gene_hits) == 0) {
    warning("plot_de_heatmap: no genes passed filtering in any comparison; returning NULL.")
    return(NULL)
  }
  
  provenance <- aggregate(
    source ~ Gene,
    data = gene_hits,
    FUN = function(x) paste(unique(x), collapse = ", ")
  )
  
  genes_wanted <- unique(gene_hits$Gene)
  
  ## compute avg_expr internally if not supplied

  if (is.null(avg_expr)) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("plot_de_heatmap: computing avg_expr internally requires the Seurat package.")
    }
    
    if (!group_by %in% colnames(seurat_obj@meta.data)) {
      stop(
        "plot_de_heatmap: `group_by` = \"", group_by, "\" not found in ",
        "colnames(seurat_obj@meta.data). Available columns: ",
        paste(colnames(seurat_obj@meta.data), collapse = ", ")
      )
    }
    
    genes_in_obj <- intersect(genes_wanted, rownames(seurat_obj))
    if (length(genes_in_obj) == 0) {
      stop("plot_de_heatmap: none of the DE genes were found in `seurat_obj`.")
    }
    
    use_assay <- if (!is.null(assay)) assay else Seurat::DefaultAssay(seurat_obj)
    
    seurat_version <- utils::packageVersion("Seurat")
    mat <- if (seurat_version >= "5.0.0") {
      SeuratObject::LayerData(seurat_obj, assay = use_assay, layer = slot)
    } else {
      Seurat::GetAssayData(seurat_obj, assay = use_assay, slot = slot)
    }
    
    mat <- mat[genes_in_obj, , drop = FALSE]
    
    groups <- seurat_obj@meta.data[[group_by]]
    if (ncol(mat) != length(groups)) {
      stop(
        "plot_de_heatmap: number of cells in the expression matrix (", ncol(mat),
        ") does not match length of `group_by` metadata (", length(groups), ")."
      )
    }
    
    group_idx <- split(seq_along(groups), groups)
    avg_expr <- sapply(group_idx, function(idx) Matrix::rowMeans(mat[, idx, drop = FALSE]))
    avg_expr <- as.data.frame(avg_expr)
  }
  
  avg_expr <- as.data.frame(avg_expr)
  
  missing_genes <- setdiff(genes_wanted, rownames(avg_expr))
  if (length(missing_genes) > 0) {
    warning(
      "plot_de_heatmap: ", length(missing_genes),
      " gene(s) from de_results not found in avg_expr and will be dropped: ",
      paste(head(missing_genes, 10), collapse = ", "),
      if (length(missing_genes) > 10) ", ..." else ""
    )
  }
  
  genes_present <- intersect(genes_wanted, rownames(avg_expr))
  
  if (length(genes_present) < 2) {
    warning("plot_de_heatmap: fewer than 2 genes available to plot after matching; returning NULL.")
    return(NULL)
  }
  
  mat <- as.matrix(avg_expr[genes_present, timepoint_order, drop = FALSE])
  
  ## score each gene across timepoints 

  mat_scaled <- t(scale(t(mat)))
  
  zero_var_genes <- rownames(mat_scaled)[apply(mat_scaled, 1, function(x) any(is.nan(x)))]
  if (length(zero_var_genes) > 0) {
    mat_scaled[zero_var_genes, ] <- 0
    warning(
      "plot_de_heatmap: ", length(zero_var_genes),
      " gene(s) had ~zero variance across timepoints in avg_expr (set to 0 in heatmap): ",
      paste(zero_var_genes, collapse = ", ")
    )
  }
  
  ## order rows by hierarchical clustering of temporal pattern
  if (nrow(mat_scaled) >= 3) {
    hc <- stats::hclust(stats::dist(mat_scaled), method = "complete")
    gene_order <- rownames(mat_scaled)[hc$order]
  } else {
    gene_order <- rownames(mat_scaled)
  }
  
  df_long <- data.frame(
    Gene = rep(rownames(mat_scaled), times = ncol(mat_scaled)),
    Timepoint = rep(colnames(mat_scaled), each = nrow(mat_scaled)),
    z = as.vector(mat_scaled),
    stringsAsFactors = FALSE
  )
  
  df_long$Gene <- factor(df_long$Gene, levels = gene_order)
  df_long$Timepoint <- factor(df_long$Timepoint, levels = timepoint_order)
  
  p <- ggplot(
    df_long,
    aes(x = Timepoint, y = Gene, fill = z)
  ) +
    
    geom_tile(colour = "white", linewidth = 0.3) +
    
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      midpoint = 0,
      name = "z-score\n(expression)"
    ) +
    
    labs(
      title = title,
      subtitle = sprintf(
        "Union of top %d up/down genes per comparison\n
        (p_val_adj < %s, |log2FC| > %s, min pct %s%%) -- %d genes shown",
        n_genes_per_comparison, fdr_cutoff, logfc_cutoff, min_pct * 100, length(gene_order)
      ),
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30"),
      axis.text.x = element_text(colour = "black", face = "bold"),
      axis.text.y = element_text(colour = "black", size = 8),
      panel.grid = element_blank()
    )
  
  attr(p, "gene_provenance") <- provenance
  
  p
}

# ==================================================
# Cell Composition Stacked Bar Chart (Broad)
# ==================================================

plot_cell_composition_broad <- function(
    plot_data
) { 
  
  p <- ggplot(
    plot_data,
    aes(time_point, fill = cell_type)
  ) +
    
    geom_bar(position = "fill") +
    scale_y_continuous(
      labels = percent_format()
    ) +
    
    labs(
      x = "Developmental stage",
      y = "Cell Proportion",
      fill = "Cell subtype",
      title = "Cell Composition",
      subtitle = paste0(settings$experiment, "_", settings$broad_celltype)
    ) +
    theme_classic()
  
  }

# =========================================================
# Cell Composition Stacked Bar Chart (Whole)
# =========================================================

plot_cell_composition_whole <- function(
    plot_data
) { 
  
  p <- ggplot(
    plot_data,
    aes(time_point, fill = broad_celltype)
  ) +
    
    geom_bar(position = "fill") +
    scale_y_continuous(
      labels = percent_format()
    ) +
    
    labs(
      x = "Developmental stage",
      y = "Cell Proportion",
      fill = "Cell subtype",
      title = "Cell Composition",
      subtitle = paste0(settings$experiment)
    ) +
    theme_classic()
  
}

