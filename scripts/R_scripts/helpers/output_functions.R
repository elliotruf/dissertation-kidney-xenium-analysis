
# =================================
# Make Output Directories
# =================================
# Create the directory structure
# for saving analysis output

make_output_dirs <- function(
    analysis,
    experiment_name,
    subdirectories = NULL
) {
  
  root_dir <- project_path(
    "results",
    analysis,
    experiment_name
  )
  
  dirs <- list(
    root = root_dir
  )
  
  if (!is.null(subdirectories)) {
    
    for (directory in subdirectories) {
      
      dirs[[directory]] <- file.path(
        root_dir,
        directory
      )
      
    }
    
  }
  
  purrr::walk(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  invisible(dirs)
  
}

# ====================================
# Save DE Figures (edgeR)
# ====================================

save_de_figures_edge_r <- function(
    de_results,
    output_dir,
    experiment_name,
    fdr_cutoff,
    logfc_cutoff,
    n_labels
) {
  
  for (name in names(de_results)) {
    
    comparison <- gsub(
      "_vs_",
      " vs ",
      name
    )
    
    title <- paste(
      experiment_name,
      comparison,
      sep = " — "
    )
    
    p <- plot_volcano_edge_r(
      de_table = de_results[[name]],
      title = title,
      fdr_cutoff = fdr_cutoff,
      logfc_cutoff = logfc_cutoff,
      n_labels = n_labels
    )
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          experiment_name,
          "_",
          name,
          "_volcano.pdf"
        )
      ),
      plot = p,
      width = 6,
      height = 5
    )
    
  }
  
}

# ===========================
# Save DE tables (edgeR)
# ===========================

save_de_tables_edge_r <- function(
    de_results,
    output_dir,
    experiment_name,
    fdr_cutoff,
    logfc_cutoff
) {
  
  for (name in names(de_results)) {
    
    file_prefix <- paste(
      experiment_name,
      name,
      sep = "_"
    )
    
    # Save complete DE results.
    write.csv(
      de_results[[name]],
      file = file.path(
        output_dir,
        paste0(
          file_prefix,
          "_DE.csv"
        )
      )
    )
    
    # Save significant DE results.
    sig <- subset(
      de_results[[name]],
      FDR < fdr_cutoff &
        abs(logFC) > logfc_cutoff
    )
    
    write.csv(
      sig,
      file = file.path(
        output_dir,
        paste0(
          file_prefix,
          "_DE_sig.csv"
        )
      )
    )
    
  }
  
}

# ====================================
# Save GO Figures (edgeR)
# ====================================

save_go_figures_edge_r <- function(
    go_results,
    output_dir,
    experiment_name,
    n_terms
) {
  
  for (name in names(go_results)) {
    
    if (is.null(go_results[[name]])) {
      next
    }
    
    comparison <- sub(
      "_(up|down)$",
      "",
      name
    )
    
    comparison <- gsub(
      "_vs_",
      " vs ",
      comparison
    )
    
    direction <- ifelse(
      grepl("_up$", name),
      "Up-regulated genes",
      "Down-regulated genes"
    )
    
    title <- paste(
      experiment_name,
      comparison,
      direction,
      sep = " — "
    )
    
    p <- plot_go_edgeR(
      go_results[[name]],
      title = title,
      n_terms = n_terms
    )
    
    if (is.null(p)) {
      next
    }
    
    n_actual <- min(
      n_terms,
      nrow(as.data.frame(go_results[[name]]))
    )
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          experiment_name,
          "_",
          name,
          "_GO.pdf"
        )
      ),
      plot = p,
      width = 10,
      height = max(
        3.5,
        0.4 * n_actual + 2
      )
    )
    
  }
  
}

# =========================
# Save GO Tables
# =========================

save_go_tables <- function(
    go_results,
    output_dir,
    experiment_name,
    condition_value = NULL
) {
  
  for (name in names(go_results)) {
    
    if (is.null(go_results[[name]])) {
      next
    }
    
    go_df <- as.data.frame(
      go_results[[name]]
    )
    
    go_df$GeneRatio <- as.character(
      go_df$GeneRatio
    )
    
    go_df$BgRatio <- as.character(
      go_df$BgRatio
    )
    
    file_prefix <- if (is.null(condition_value)) {
      
      paste(
        experiment_name,
        name,
        sep = "_"
      )
      
    } else {
      
      paste(
        experiment_name,
        condition_value,
        name,
        sep = "_"
      )
      
    }
    
    writexl::write_xlsx(
      go_df,
      path = file.path(
        output_dir,
        paste0(
          file_prefix,
          "_GO.xlsx"
        )
      )
    )
  }
}

# ======================================
# Save DE Figures (Seurat)
# ======================================
save_de_figures_seurat <- function(
    de_results,
    output_dir,
    experiment_name,
    condition_value,
    fdr_cutoff,
    logfc_cutoff
) {
  for (name in names(de_results)) {
    comparison <- gsub("_vs_", " vs ", name)
    
    title <- paste(
      experiment_name,
      comparison,
      sep = " — "
    )
    
    p <- plot_top_genes_bar(
      de_table = de_results[[name]],
      title = title,
      fdr_cutoff = fdr_cutoff,
      logfc_cutoff = logfc_cutoff
    )
    
    if (is.null(p)) next
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          experiment_name,
          "_",
          name,
          "_DEG_bar_plot.pdf"
        )
      ),
      plot = p,
      width = 7,
      height = 5
    )
  }
}

# ======================================
# Save GO Figures (Seurat)
# ======================================
save_go_figures_seurat <- function(
    go_results,
    output_dir,
    experiment_name,
    condition_value,
    n_terms = 15
) {
  for (name in names(go_results)) {
    if (is.null(go_results[[name]])) next
    
    comparison <- gsub("_vs_", " vs ", name)
    
    direction <- ifelse(
      grepl("_up$", name),
      "Up-regulated genes",
      "Down-regulated genes"
    )
    
    title <- paste(
      experiment_name,
      comparison,
      direction,
      sep = " — "
    )
    
    p <- plot_go_seurat(
      go_result = go_results[[name]],
      title = title,
      n_terms = n_terms
    )
    
    if (is.null(p)) next
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          experiment_name,
          "_",
          name,
          "_GO.pdf"
        )
      ),
      plot = p,
      width = 10,
      height = 7
    )
  }
}

# ==========================
# Save DE Tables (Seurat)
# ==========================
save_de_tables_seurat <- function(
    de_results,
    output_dir,
    experiment_name,
    condition_value,
    fdr_cutoff,
    logfc_cutoff
) {
  for (name in names(de_results)) {
    
    file_prefix <- paste(
      experiment_name,
      name,
      sep = "_"
    )
    
    write.csv(
      de_results[[name]],
      file = file.path(
        output_dir,
        paste0(file_prefix, "_DE.csv")
      )
    )
    
    sig <- subset(
      de_results[[name]],
      p_val_adj < fdr_cutoff &
        abs(avg_log2FC) > logfc_cutoff
    )
    
    write.csv(
      sig,
      file = file.path(
        output_dir,
        paste0(file_prefix, "_DE_sig.csv")
      )
    )
  }
}

# ======================================
# Save DE Heatmap (Seurat)
# ======================================

# ======================================
# Save DE Heatmap (Seurat)
# ======================================
save_de_heatmap_seurat <- function(
    de_results,
    timepoint_order,
    output_dir,
    experiment_name,
    cell_type,
    condition_value,
    seurat_obj,
    group_by,
    assay,
    fdr_cutoff,
    logfc_cutoff,
    min_pct = 0.1,
    n_genes_per_comparison = 20,
    slot = "data",
    width = 7,
    height = 10
) {
  title <- paste(
    experiment_name,
    cell_type,
    condition_value,
    "DE Genes Across Timecourse",
    sep = " — "
  )
  
  p <- plot_de_heatmap(
    de_results = de_results,
    seurat_obj = seurat_obj,
    group_by = group_by,
    assay = assay,
    slot = slot,
    timepoint_order = timepoint_order,
    title = title,
    fdr_cutoff = fdr_cutoff,
    logfc_cutoff = logfc_cutoff,
    min_pct = min_pct,
    n_genes_per_comparison = n_genes_per_comparison
  )
  
  if (is.null(p)) {
    warning(
      "save_de_heatmap_seurat: plot_de_heatmap() returned NULL; ",
      "heatmap not saved."
    )
    return(invisible(NULL))
  }
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        experiment_name,
        "_",
        cell_type,
        "_",
        condition_value,
        "_DEG_heatmap.pdf"
      )
    ),
    plot = p,
    width = width,
    height = height
  )
  
  invisible(p)
}

# =================================
# Save Cell Composition Tables
# =================================

save_composition_summary_tables <- function(
    summary_list, 
    output_dir
) { 
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  purrr::iwalk(
    
    summary_list,
    
    \(tbl, nm)
    
    readr::write_csv(
      tbl,
      file.path(output_dir, paste0(nm, ".csv"))
    )
  )
  
  }
