# Script for the extraction of ROIs and generation of Seurat objects.

# ==================
# User Parameters
# ==================

settings <- list(
  # Specify directory of ROIs
  roi_dir = "D:/USERS/ELLIOT/Dissertation/xen1/Xen1_Male_Vessels",
  
  experiment = "Xen1_Male_Vessels",
  
  # Specify order of metadata in the roi folder name
  ## Timepoints handled by child directories, not here.
  ## Group refers to sex (xen1), genotype (xen2), etc.
  metadata = c("sample", "group", "tissue"),
  
  seurat_obj = xen1diet,
  
  timepoint_pattern = "(?i)naive|sham|\\d+h|\\d+d|wk\\d+|\\d+wk"
)   

# ==================

# Modules
library(Seurat)
library(readxl)
library(tidyverse)
library(fs)

source("scripts/helpers/extraction_functions_v2.R")

# Get ROIs from ROI directory
roi_files <- find_rois(settings$roi_dir)

# Build ROI tibble
roi_tibble <- build_roi_tibble(roi_files)

# Build ROI Seurat Objects
roi_objects <- build_roi_objects(roi_tibble, settings$seurat_obj)

# Save ROI objects
saveRDS(
  roi_objects,
  file = file.path(
    "results",
    "objects",
    paste0(settings$experiment, "_roi_objects.rds")
)
)

# Save ROI tibble
saveRDS(
  roi_tibble,
  file = file.path(
    "results",
    "objects",
    paste0(settings$experiment, "_roi_tibble.rds")
  )
)
