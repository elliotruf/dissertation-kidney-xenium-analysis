# ====================================
# Function 1: Save DE figures (EdgeR)
# ====================================

save_de_figures_edgeR <- function(
    de_results,
    output_dir,
    settings,
    run_name = NULL
) {
  
  for (name in names(de_results)) {
    
    comparison <- gsub("_vs_", " vs ", name)
    
    run_label <- if (!is.null(run_name)) {
      gsub("_", " ", run_name)
    } else {
      NULL
    }
    
    title <- if (is.null(run_label)) {
      comparison
    } else {
      paste(run_label, "—", comparison)
    }
    
    p <- plot_volcano_edgeR(
      de_table = de_results[[name]],
      title = title,
      fdr_cutoff = settings$fdr_cutoff,
      logfc_cutoff = settings$logfc_cutoff,
      n_labels = settings$n_labels
    )
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(name, "_volcano.pdf")
      ),
      plot = p,
      width = 6,
      height = 5
    )
  }
}

# ====================================
# Function 2: Save GO figures (EdgeR)
# ====================================

save_go_figures_edgeR <- function(
    go_results,
    output_dir,
    settings,
    run_name = NULL
) {
  
  for (name in names(go_results)) {
    
    comparison <- sub("_(up|down)$", "", name)
    comparison <- gsub("_vs_", " vs ", comparison)
    
    direction <- ifelse(
      grepl("_up$", name),
      "Up-regulated genes",
      "Down-regulated genes"
    )
    
    run_label <- if (!is.null(run_name)) {
      gsub("_", " ", run_name)
    } else {
      NULL
    }
    
    title <- if (is.null(run_label)) {
      paste(comparison, direction, sep = " — ")
    } else {
      paste(run_label, comparison, direction, sep = " — ")
    }
    
    p <- plot_go_edgeR(
      go_results[[name]],
      title = title,
      n_terms = settings$go_n_terms
    )
    
    if (is.null(p))
      next
    
    ggsave(
      
      filename = file.path(
        output_dir,
        paste0(
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

# ===========================
# Function 3: Save DE tables (EdgeR)
# ===========================

save_de_tables_edgeR <- function(
    de_results,
    output_dir,
    settings
) {
  
  for (name in names(de_results)) {
    
    write.csv(
      de_results[[name]],
      
      file.path(
        output_dir,
        paste0(name, "_DE.csv")
      )
      
    )
    
    sig <- subset(
      
      de_results[[name]],
      
      FDR < settings$fdr_cutoff &
        abs(logFC) > settings$logfc_cutoff
      
    )
    
    write.csv(
      sig,
      
      file.path(
        output_dir,
        paste0(name, "_DE_sig.csv")
      )
      
    )
    
  }
  
}

# =========================
# Function 4: Save GO tables
# =========================

save_go_tables <- function(
    go_results,
    output_dir
) {
  
  for (name in names(go_results)) {
    
    if (is.null(go_results[[name]]))
      next
    
    write.csv(
      as.data.frame(go_results[[name]]),
      
      file.path(
        output_dir,
        paste0(name, "_GO.csv")
      ),
      
      row.names = FALSE
      
    )
    
  }
  
}

# =================================
# Function 5: Output Directories
# =================================
# Create the directory structure
# for saving analysis output

make_output_dirs <- function(
    analysis,
    run_name
) {
  
  dirs <- list(
    
    de_tables = project_path(
      "results",
      analysis,
      "DE",
      "tables",
      run_name
    ),
    
    de_figures = project_path(
      "results",
      analysis,
      "DE",
      "figures",
      run_name
    ),
    
    go_tables = project_path(
      "results",
      analysis,
      "GO",
      "tables",
      run_name
    ),
    
    go_figures = project_path(
      "results",
      analysis,
      "GO",
      "figures",
      run_name
    )
    
  )
  
  purrr::walk(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
    
  )
  
  invisible(dirs)
  
}

# ======================================
# Function 6: Save DE Figures (Seurat)
# ======================================

save_de_figures_seurat <- function(
    de_results,
    output_dir,
    settings,
    run_name = NULL
) {
  
  for (name in names(de_results)) {
    
    comparison <- gsub("_vs_", " vs ", name)
    
    run_label <- if (!is.null(run_name)) {
      gsub("_", " ", run_name)
    } else {
      NULL
    }
    
    title <- if (is.null(run_label)) {
      comparison
    } else {
      paste(run_label, "—", comparison)
    }
    
    p <- plot_top_genes_bar(
      de_table = de_results[[name]],
      title    = name,
      fdr_cutoff = settings$fdr_cutoff,
      logfc_cutoff = settings$logfc_cutoff
    )

    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(name, "_DEG_bar_plot.pdf")
      ),
      plot = p,
      width = 7,
      height = 5
    )
    
  }
  
}


# ======================================
# Function 7: Save GO Figures (Seurat)
# ======================================

save_go_figures_seurat <- function(
    go_results,
    output_dir,
    settings,
    run_name = NULL
) {
  
  for (name in names(go_results)) {
    
    comparison <- sub("_(up|down)$", "", name)
    comparison <- gsub("_vs_", " vs ", comparison)
    
    direction <- ifelse(
      grepl("_up$", name),
      "Up-regulated genes",
      "Down-regulated genes"
    )
    
    run_label <- if (!is.null(run_name)) {
      gsub("_", " ", run_name)
    } else {
      NULL
    }
    
    title <- if (is.null(run_label)) {
      paste(comparison, direction, sep = " — ")
    } else {
      paste(run_label, comparison, direction, sep = " — ")
    }
    
    p <- plot_go_seurat(
      go_results[[name]],
      title = title,
      n_terms = settings$go_n_terms
    )
    
    if (is.null(p))
      next
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(name, "_GO.pdf")
      ),
      plot = p,
      width = 10,
      height = 7
    )
    
  }
  
}

# ====================================
# Function 8: Save DE Tables (Seurat)
# ====================================

save_de_tables_seurat <- function(
    de_results,
    output_dir,
    settings
) {

  for (name in names(de_results)) {

    write.csv(
      de_results[[name]],
      file.path(
        output_dir,
        paste0(name, "_DE.csv")
      )
    )

    sig <- subset(

      de_results[[name]],

      p_val_adj < settings$fdr_cutoff &
        abs(avg_log2FC) > settings$logfc_cutoff

    )

    write.csv(
      sig,
      file.path(
        output_dir,
        paste0(name, "_DE_sig.csv")
      )
    )

  }

}

# ======================================
# Function 9: Save DE Heatmap (Seurat)
# ======================================

save_de_heatmap_seurat <- function(
    de_results,
    timepoint_order,
    output_dir,
    settings,
    run_name = NULL,
    avg_expr = NULL,     # supply this OR seurat_obj/group_by (see plot_de_heatmap)
    seurat_obj = NULL,
    group_by = NULL,
    assay = NULL,
    slot = "data",
    width = 7,
    height = 10
) {
  
  run_label <- if (!is.null(run_name)) gsub("_", " ", run_name) else NULL
  
  title <- if (is.null(run_label)) {
    "DE Genes Across Timecourse"
  } else {
    paste(run_label, "\u2014 DE Genes Across Timecourse")
  }
  
  p <- plot_de_heatmap(
    de_results   = de_results,
    avg_expr     = avg_expr,
    seurat_obj   = seurat_obj,
    group_by     = group_by,
    assay        = assay,
    slot         = slot,
    timepoint_order = timepoint_order,
    title        = title,
    fdr_cutoff   = settings$fdr_cutoff,
    logfc_cutoff = settings$logfc_cutoff,
    min_pct      = if (!is.null(settings$min_pct)) settings$min_pct else 0.1,
    n_genes_per_comparison = if (!is.null(settings$n_genes_per_comparison)) {
      settings$n_genes_per_comparison
    } else {
      20
    }
  )
  
  if (is.null(p)) {
    warning(
      "save_de_heatmap_seurat: plot_de_heatmap() returned NULL ",
      "(no genes passed filtering across any comparison); heatmap not saved."
    )
    return(invisible(NULL))
  }
  
  out_name <- paste0(
    if (!is.null(run_name)) paste0(run_name, "_") else "",
    "DEG_heatmap.pdf"
  )
  
  ggsave(
    filename = file.path(output_dir, out_name),
    plot     = p,
    width    = width,
    height   = height
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
