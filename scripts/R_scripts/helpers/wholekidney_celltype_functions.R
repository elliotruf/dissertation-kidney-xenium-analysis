
# =================================
# Define Broad Celltypes
# =================================
# Collapse fine cell type annotations into
# broader types for downstream analyses

assign_broad_celltypes <- function(cell_type) {
  
  case_when(
    grepl("^PTS", cell_type) ~ "PT",
    grepl("^TAL", cell_type) ~ "TAL",
    grepl("^ATL|^DTL", cell_type) ~ "Thin limb",
    grepl("^CD", cell_type) ~ "Collecting duct",
    grepl("^DCT", cell_type) ~ "DCT",
    grepl("^STROMA", cell_type) ~ "Stroma",
    grepl("^ENDO", cell_type) ~ "Endothelium",
    grepl("^GEnC", cell_type) ~ "Endothelium",
    grepl("^VSMC", cell_type) ~ "Vascular",
    grepl("^IMMUNE", cell_type) ~ "Immune",
    grepl("^MESANGIAL", cell_type) ~ "Glomerular",
    grepl("^PODO", cell_type) ~ "Glomerular",
    grepl("^PEC", cell_type) ~ "Glomerular",
    grepl("^UROTHELIUM", cell_type) ~ "Urothelium",
    grepl("^FAT", cell_type) ~ "Adipose",
    grepl("^rare", cell_type) ~ "Rare",
    TRUE ~ "Other"
  )
}

# ==========================================
# Function 2: Add Broad Cell Types
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
# Function 3: Subset Broad Cell Type
# ==========================================
# Extract a single broad cell type for
# downstream analysis.

subset_celltype <- function(
    seurat_obj,
    celltype,
    column = "broad_celltype"
) {
  
  if (!column %in% colnames(seurat_obj@meta.data))
    stop(paste("Column", column, "not found."))
  
  cells <- rownames(
    seurat_obj@meta.data[
      seurat_obj@meta.data[[column]] == celltype,
      ,
      drop = FALSE
    ]
  )
  
  subset(
    seurat_obj,
    cells = cells
  )
  
}

# ==========================================
# Function 4: Cell Count Summary
# ==========================================
# Summarise numbers of cells by condition.

summarise_celltype <- function(
    seurat_obj,
    sample_col = "sample_id",
    variable = "time_point"
) {
  
  seurat_obj@meta.data %>%
    count(
      .data[[sample_col]],
      .data[[variable]],
      name = "cells"
    ) %>%
    group_by(.data[[sample_col]]) %>%
    mutate(
      proportion = cells / sum(cells)
    ) %>%
    ungroup()
  
}

# =================================================================
# Function 5: Whole-kidney Broad Cell-type Differential Expression 
# =================================================================
# Compare two experimental groups using
# Seurat FindMarkers.

run_celltype_de <- function(
    seurat_obj,
    ident1,
    ident2,
    variable,
    assay = "Xenium"
) {
  
  DefaultAssay(seurat_obj) <- assay
  
  Idents(seurat_obj) <- seurat_obj[[variable, drop = TRUE]]
  
  FindMarkers(
    seurat_obj,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
}

# ==========================================
# Function 7: Run GO For Seurat DE Results
# ==========================================
# Perform GO enrichment separately for
# up- and down-regulated genes from
# Seurat FindMarkers() output.

run_all_go_seurat <- function(
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
      subset(sig, avg_log2FC > 0)
    )
    
    down_genes <- rownames(
      subset(sig, avg_log2FC < 0)
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

# ==========================================
# Function 8: Run Name
# ==========================================
# Generate a unique identifier for the
# current cell type analysis.

get_celltype_run_name <- function(
    settings
) {
  
  paste(
    settings$experiment,
    settings$celltype,
    sep = "_"
  )
  
}