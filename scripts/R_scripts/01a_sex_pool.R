# Script for combining ("pooling") two .rds objects
# such as to combine Male and Female ROI selections for analysis
# Particularly useful for ROI selections of minimal cells

# ===================
# User Settings
# ===================

settings <- list(
  
  seurat_obj_1 = file.path(
    "results",
    "objects",
    "Xen1_Male_Vessels_roi_objects.rds"
  ),
  
  seurat_obj_2 = file.path(
    "results",
    "objects",
    "Xen1_Female_Vessels_roi_objects.rds"
  ),
  
  experiment = "Xen1",
  
  tissue = "Vessels"
    
)

# ===================
# Dependencies
# ===================

library(Seurat)

source("scripts/R_scripts/helpers/project_paths_v3.R")

# ================
# Pool Objects
# ================

seurat_obj_1 <- readRDS(project_path(settings$seurat_obj_1))

seurat_obj_2 <- readRDS(project_path(settings$seurat_obj_2))

message("Pooling objects...")
seurat_obj_combined <- c(
  seurat_obj_1,
  seurat_obj_2
)

# ========================
# Save Combined Object
# ========================

message("Saving Combined Object...")
seurat_obj_combined <- saveRDS(
  seurat_obj_combined,
  
  file = project_path(
    "results",
    "objects",
    
    paste0(
      settings$experiment, "_",
      settings$tissue, "_",
      "combined_roi_objects.rds"
    )
  )
)
  

