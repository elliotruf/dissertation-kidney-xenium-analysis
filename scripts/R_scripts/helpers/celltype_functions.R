
# =============================================
# Function 1: Define Broad Celltypes (Custom)
# =============================================
# Collapse fine cell type annotations into
# broader types for downstream analyses

assign_broad_celltypes <- function(cell_type) {
  
  dplyr::case_when(
    
    # Proximal tubule
    grepl("^PTS", cell_type) |
      grepl("^PT_", cell_type) |
      cell_type %in% c(
        "Prolif_PT",
        "NewPT2_PEC"
      ) ~ "PT",
    
    # TAL
    grepl("^TAL", cell_type) |
      cell_type %in% c("Prolif_TAL") ~ "TAL",
    
    # Thin limb
    grepl("^ATL", cell_type) |
      grepl("^DTL", cell_type) ~ "Thin limb",
    
    # DCT
    grepl("^DCT", cell_type) |
      cell_type %in% c(
        "CNT",
        "CNT _INJURY",
        "DCT_CNT",
        "Prolif_DT"
      ) ~ "DCT",
    
    # Collecting duct
    grepl("^CD", cell_type) |
      cell_type %in% c(
        "IC_A",
        "IC_B",
        "PC_OMCD",
        "PC_IMCD"
      ) ~ "Collecting duct",
    
    # Endothelium
    grepl("^ENDO", cell_type) |
      grepl("^GEnC", cell_type) |
      cell_type %in% c(
        "EC",
        "Prolif_EC",
        "Endothelial_vessels"
      ) ~ "Endothelium",
    
    # Stroma
    grepl("^STROMA", cell_type) |
      cell_type %in% c(
        "Fib",
        "MC",
        "Prolif_Fib"
      ) ~ "Stroma",
    
    # Vascular
    grepl("^VSMC", cell_type) |
      cell_type %in% c(
        "VSMCs",
        "VSMCs_small_vessels",
        "VSMCs_large vessels",
        "VSMCs_mature_M"
      ) ~ "Vascular",
    
    # Immune
    grepl("^IMMUNE", cell_type) |
      cell_type %in% c(
        "Immune",
        "Prolif_Immune"
      ) ~ "Immune",
    
    # Glomerular
    grepl("^MESANGIAL", cell_type) |
      grepl("^PEC", cell_type) |
      grepl("^PODO", cell_type) |
      grepl("^Podocytes", cell_type) |
      cell_type %in% c(
        "JGC",
        "Pod"
      ) ~ "Glomerular",
    
    # Urothelium
    grepl("^UROTHELIUM", cell_type) |
      grepl("^Urothelium", cell_type) |
      cell_type %in% c("Uro") ~ "Urothelium",
    
    # Adipose
    grepl("^FAT", cell_type) ~ "Adipose",
    
    # Rare
    grepl("^rare", cell_type) |
      grepl("^Rare", cell_type) ~ "Rare",
    
    # Unresolved
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
    
    seurat_obj$broad_celltype <-
      
      assign_broad_celltypes(
        seurat_obj[[fine_column, drop = TRUE]]
      )
  
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
# Function 4: Whole Dataset Summary
# ==========================================
# Summarise whole-dataset cell numbers and
# broad/fine cell-type composition by condition.

summarise_whole_data <- function(
    seurat_obj,
    condition_col = "sex",
    time_col = "time_point",
    broad_col = "broad_celltype",
    fine_col = "cell_type"
) {
  
  meta <- seurat_obj[[]]
  
  # Check required columns exist
  required_cols <- c(
    condition_col,
    time_col,
    broad_col,
    fine_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(meta)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing metadata columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Remove cells with missing condition or time-point metadata
  meta <- meta %>%
    filter(
      !is.na(.data[[condition_col]]),
      !is.na(.data[[time_col]])
    )
  
  # Total cells per condition and time point
  dataset_summary <- meta %>%
    group_by(
      .data[[condition_col]],
      .data[[time_col]]
    ) %>%
    summarise(
      total_cells = n(),
      .groups = "drop"
    )
  
  # Broad cell-type composition
  whole_composition_broad <- meta %>%
    count(
      .data[[condition_col]],
      .data[[time_col]],
      .data[[broad_col]],
      name = "cells"
    ) %>%
    group_by(
      .data[[condition_col]],
      .data[[time_col]]
    ) %>%
    mutate(
      proportion = cells / sum(cells)
    ) %>%
    ungroup()
  
  # Fine cell-type composition
  whole_composition_fine <- meta %>%
    count(
      .data[[condition_col]],
      .data[[time_col]],
      .data[[fine_col]],
      name = "cells"
    ) %>%
    group_by(
      .data[[condition_col]],
      .data[[time_col]]
    ) %>%
    mutate(
      proportion = cells / sum(cells)
    ) %>%
    ungroup()
  
  # Fine-to-broad annotation lookup
  lookup <- meta %>%
    distinct(
      .data[[fine_col]],
      .data[[broad_col]]
    )
  
  list(
    dataset_summary = dataset_summary,
    whole_composition_broad = whole_composition_broad,
    whole_composition_fine = whole_composition_fine,
    celltype_lookup = lookup
  )
}

# ==========================================
# Function: ROI Summary
# ==========================================
# Summarise cell numbers and composition
# within ROI-specific Seurat objects.

summarise_rois <- function(
    roi_objects,
    roi_col = "roi",
    condition_col = "sex",
    tissue_col = "tissue",
    time_col = "time_point",
    broad_col = "broad_celltype",
    fine_col = "cell_type"
) {
  
  meta <- purrr::map_dfr(
    roi_objects,
    ~ .x@meta.data
  )
  
  # Check required metadata columns
  required_cols <- c(
    roi_col,
    condition_col,
    tissue_col,
    time_col,
    broad_col,
    fine_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(meta)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing metadata columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Total cells per ROI
  roi_summary <-
    meta %>%
    group_by(
      .data[[roi_col]],
      .data[[condition_col]],
      .data[[tissue_col]],
      .data[[time_col]]
    ) %>%
    summarise(
      total_cells = n(),
      .groups = "drop"
    )
  
  # Broad cell type composition per ROI
  roi_composition_broad <-
    meta %>%
    count(
      .data[[roi_col]],
      .data[[condition_col]],
      .data[[tissue_col]],
      .data[[time_col]],
      .data[[broad_col]],
      name = "cells"
    ) %>%
    group_by(
      .data[[roi_col]]
    ) %>%
    mutate(
      proportion = cells / sum(cells)
    ) %>%
    ungroup()
  
  # Fine cell type composition per ROI
  roi_composition_fine <-
    meta %>%
    count(
      .data[[roi_col]],
      .data[[condition_col]],
      .data[[tissue_col]],
      .data[[time_col]],
      .data[[fine_col]],
      name = "cells"
    ) %>%
    group_by(
      .data[[roi_col]]
    ) %>%
    mutate(
      proportion = cells / sum(cells)
    ) %>%
    ungroup()
  
  list(
    roi_summary = roi_summary,
    roi_composition_broad = roi_composition_broad,
    roi_composition_fine = roi_composition_fine
  )
  
}

# ==========================================
# Function: Differential Expression
# ==========================================
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
# Function: Run GO For Seurat DE Results
# ==========================================
# Perform GO enrichment separately for
# up- and down-regulated genes from
# Seurat FindMarkers() output.

run_all_go_seurat <- function(
    de_results,
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
# Function: Run Name
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

# ===================================================
# Function: Assign Broad Cell Types (XeniumClean)
# ===================================================

assign_broad_celltypes_xeniumclean <- function(cell_type) {
  
  cell_type <- as.character(cell_type)
  
  dplyr::case_when(
    cell_type %in% c(
      "PTS1", "PTS1S2", "PTS2", "PTS2S3", "PTS3",
      "PT_Immature", "NewPT2_PEC"
    ) ~ "PT",
    
    cell_type %in% c("DTL", "ATL") ~ "Thin limb",
    
    cell_type %in% c("TAL", "TAL_DCT", "Prolif_TAL") ~ "TAL",
    
    cell_type %in% c("DCT1", "DCT2", "CNT") ~ "DCT",
    
    cell_type %in% c(
      "PC_IMCD", "PC_OMCD",
      "IC_A", "IC_B",
      "Prolif_DT"
    ) ~ "Collecting duct",
    
    cell_type %in% c("Pod", "JGC") ~ "Glomerular",
    
    cell_type %in% c("EC", "Prolif_EC") ~ "Endothelium",
    
    cell_type %in% c("Fib", "MC", "Prolif_Fib") ~ "Stroma",
    
    cell_type %in% c("Immune", "Prolif_Immune") ~ "Immune",
    
    cell_type == "Uro" ~ "Urothelium",
    
    TRUE ~ "Other"
  )
}

