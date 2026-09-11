# Script for combining ("pooling") two .rds objects
# such as to combine Male and Female ROI selections for analysis
# Particularly useful for ROI selections of minimal cells

# Sex Pooling
#
# Pools two Seurat objects together for combined analysis
# Particularly useful for ROI selections with minimal cells
# Such as male and female vessels
#
# Before running:
# - Set the Seurat objects to be merged.
# - Set the experiment name.
# - Set the tissue name

# ===================
# User Settings
# ===================

settings <- list(
  
  seurat_obj_1 = file.path(
    "results",
    "ROI_extraction",
    "Xen1_Male_Vessels",
    "objects",
    "Xen1_Male_Vessels_roi_objects.rds"
  ),
  
  seurat_obj_2 = file.path(
    "results",
    "ROI_extraction",
    "Xen1_Female_Vessels",
    "objects",
    "Xen1_Female_Vessels_roi_objects.rds"
  ),
  
  experiment_name = "Xen1",
  
  tissue = "Vessels"
    
)

# ===================
# Dependencies
# ===================

library(Seurat)

source("scripts/R_scripts/helpers/project_paths.R")

# ================
# Pool Objects
# ================

message("Loading objects...")
seurat_obj_1 <- readRDS(project_path(settings$seurat_obj_1))

seurat_obj_2 <- readRDS(project_path(settings$seurat_obj_2))

message("Pooling objects...")
seurat_obj_combined <- c(
  seurat_obj_1,
  seurat_obj_2
)

# ========================
# Make Output Directory
# ========================

output_dirs <- make_output_dirs(
  analysis = "ROI_extraction",
  experiment_name = paste0(
    settings$experiment_name,
    "_",
    settings$tissue,
    "_combined"
  ),
  subdirectories = c(
    "objects"
  )
)

message("Saving results...")

# ========================
# Save Combined Object
# ========================

message("Saving Combined Object...")

saveRDS(
  seurat_obj_combined,
  file = file.path(
    output_dirs$object,
    paste0(
      settings$experiment_name,
      "_",
      settings$tissue,
      "_combined_roi_objects.rds"
    )
  )
)

# Done!
message("Done!")
