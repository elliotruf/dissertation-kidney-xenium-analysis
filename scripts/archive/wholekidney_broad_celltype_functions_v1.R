# ========================================================
# Functions for Whole Kidney Cell Type Analysis Pipeline
# ========================================================

# ===========================
# Function 1: Get Run Name
# ===========================

get_run_name <- function(seurat_obj, object_path) {
  tools::file_path_sans_ext(basename(object_path))
}

# ====================================
# Function 2: Run Broad Cell Type DE
# ====================================

run_broadtype_de <- function(
    seurat_obj,
    broad_type,
    settings
) {
  
  obj <- subset(
    seurat_obj,
    subset = broad_type == !!broad_type
  )
  
  Idents(obj) <- obj[[settings$variable]][,1]
  
  groups <- setdiff(
    levels(Idents(obj)),
    settings$reference
  )
  
  results <- purrr::map(
    groups,
    function(g) {
      
      message("Running comparison: ", g)
      
      de <- FindMarkers(
        obj,
        ident.1 = g,
        ident.2 = settings$reference,
        assay = settings$assay,
        slot = "data"
      )
      
      ## overwrite Seurat's fold change with a manually calculated one
      cells1 <- WhichCells(obj, idents = g)
      cells2 <- WhichCells(obj, idents = settings$reference)
      
      de$avg_log2FC <- calc_logFC(
        obj,
        cells1 = cells1,
        cells2 = cells2,
        assay = settings$assay,
        layer = "data"
      )[rownames(de)]
      
      de
    }
  )
  
  names(results) <-
    paste(
      broad_type,
      groups,
      "vs",
      settings$reference,
      sep = "_"
    )
  
  results
}

# =====================
# Run GO Enrichment
# =====================

run_go_seurat <- function(
  de_results,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
) {
  
  go_results <- list()
  
  for (name in names(de_results)) {
    
    de <- de_results[[name]]
    
    # Upregulated genes
    up <- rownames(
      subset(
        de,
        p_val_adj < p_cutoff &
          avg_log2FC > 0
      )
    )
    
    if (length(up) >= min_genes) {
      
      go_results[[paste0(name, "_up")]] <-
        run_go(
          up,
          ontology = ontology
        )
      
    }
    
    # Downregulated genes
    down <- rownames(
      subset(
        de,
        p_val_adj < p_cutoff &
          avg_log2FC < 0
      )
    )
    
    if (length(down) >= min_genes) {
      
      go_results[[paste0(name, "_down")]] <-
        run_go(
          down,
          ontology = ontology
        )
      
    }
    
  }
  
  go_results
}

# =======================
# Function: Calculate LogFC
# =======================

calc_logFC <- function(
    object,
    cells1,
    cells2,
    assay = "Xenium",
    layer = "data",
    pseudocount = 1e-6
){
  
  mat <- GetAssayData(
    object,
    assay = assay,
    layer = layer
  )
  
  mean1 <- Matrix::rowMeans(mat[, cells1, drop = FALSE])
  mean2 <- Matrix::rowMeans(mat[, cells2, drop = FALSE])
  
  log2((mean1 + pseudocount)/(mean2 + pseudocount))
  
}
