# ====================================================
# Functions for the whole-kidney cell type analysis
# ====================================================

# ==========================================
# Add Broad Cell Types
# ==========================================
# Add broad cell type annotations if they
# are not already present.

add_broad_celltypes <- function(
    seurat_obj,
    fine_column = "cell_type"
) {
  
  if (!"broad_celltype" %in%
      colnames(seurat_obj@meta.data)) {
    
    seurat_obj$broad_celltype <-
      
      assign_broad_celltypes(
        seurat_obj[[fine_column, drop = TRUE]]
      )
    
  }
  
  seurat_obj
  
}


# ==========================================
# Run Cell-type Differential Expression
# ==========================================
# Compare two groups using Seurat FindMarkers.

run_cell_type_de <- function(
    seurat_obj,
    group_1,
    group_2,
    group_variable,
    assay = "Xenium"
) {
  
  DefaultAssay(seurat_obj) <- assay
  
  Idents(seurat_obj) <-
    seurat_obj[[group_variable, drop = TRUE]]
  
  FindMarkers(
    seurat_obj,
    ident.1 = group_1,
    ident.2 = group_2,
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
}

# ==========================================
# Run GO for All Seurat DE Results
# ==========================================
# Perform GO enrichment separately for
# up- and down-regulated genes from
# Seurat FindMarkers() results.

run_go_all_seurat <- function(
    de_results,
    universe,
    fdr_cutoff = 0.05,
    ontology = "BP",
    p_cutoff = 0.05,
    min_genes = 10
) {
  
  go_results <- list()
  
  for (comparison in names(de_results)) {
    
    de <- de_results[[comparison]]
    
    sig <- subset(
      de,
      p_val_adj < fdr_cutoff
    )
    
    up_genes <- rownames(
      subset(
        sig,
        avg_log2FC > 0
      )
    )
    
    down_genes <- rownames(
      subset(
        sig,
        avg_log2FC < 0
      )
    )
    
    go_results[[paste0(comparison, "_up")]] <-
      if (length(up_genes) >= min_genes) {
        
        run_go(
          genes = up_genes,
          universe = universe,
          ontology = ontology,
          p_cutoff = p_cutoff
        )
        
      } else {
        
        NULL
        
      }
    
    go_results[[paste0(comparison, "_down")]] <-
      if (length(down_genes) >= min_genes) {
        
        run_go(
          genes = down_genes,
          universe = universe,
          ontology = ontology,
          p_cutoff = p_cutoff
        )
        
      } else {
        
        NULL
        
      }
    
  }
  
  go_results
  
}
