# ROI Extraction

# ==================================
# User Settings
# ==================================

settings <- list(
  
  # Directory containing ROI folders, starting from "Dissertation"
  roi_dir = file.path(
    "data",
    "Xen2_ROIs",
    "Xen2_KO_Vessels"
  ),
  
  # Parent Seurat object
  seurat_obj = xen2diet,
  
  # Name used when saving outputs
  experiment = "Xen2_KO_Vessels",
  
  # Pattern used to extract timepoints from folder names
  time_point_pattern =
    "Naive|naive|Sham|sham|\\d+h|\\d+d|\\d+wk|wk\\d+",
  
  # Order of metadata in sample folder names
  metadata = c(
    "experiment",
    "group",
    "tissue"
  )
)

# ===============================
# Packages
# ===============================

library(Seurat)
library(tidyverse)
library(readxl)
library(fs)

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/extraction_functions_v3.R")

# ===============================
# Extract ROIs
# ===============================

message("Extracting ROIs from specified dataset...")
extract_rois(settings$seurat_obj)

# Done
message("Done!")
