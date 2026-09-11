# ROI Extraction
#
# Extracts cells belonging to each Xenium ROI from a parent Seurat object
# and saves the resulting ROI metadata and Seurat objects for downstream analysis.
#
# Before running:
# - Set the ROI directory containing the ROI "_cells_stats.csv" files.
# - set the parent Seurat object
# - Set the experiment name used for output files
# - Check the time-point pattern matches the ROI folder names.
# - Set sample_map if ROI condition labels differ from
# the parent Seurat sample_id values.

# ==================================
# User Settings
# ==================================

settings <- list(
  
  # Directory containing ROI folders
  roi_dir = file.path(
    "data",
    "Xen1_ROIs",
    "Xen1_Male_Vessels"
  ),
  
  # Parent Seurat object
  seurat_obj = file.path(
    "data",
    "xen1diet.rds"
  ),
  
  # Name used for output directories and files
  experiment_name = "Xen1_Male_Vessels",
  
  # Pattern used to extract time points from ROI folder names
  # Currently set to find time points by day, hour, or week, as well as Naive/naive or Sham/sham.
  time_point_pattern =
    "Naive|naive|Sham|sham|\\d+h|\\d+d|\\d+wk|wk\\d+",

  # Map ROI condition labels to parent Seurat sample_id values.
  # Set to NULL if labels already match.
  sample_map = NULL
  
)

# ===============================
# Dependencies
# ===============================

library(Seurat)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/ROI_extraction_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load data
# ===============================

message("Loading Seurat object...")

seurat_obj <- readRDS(settings$seurat_obj)

# ===============================
# Find ROI files
# ===============================

message("Finding ROI files from specified directory...")

roi_files <- find_roi_files(
  roi_dir = settings$roi_dir
)

# ===============================
# Build ROI metadata
# ===============================

message("Building ROI metadata...")

roi_metadata <- build_roi_metadata(
  roi_files = roi_files,
  time_point_pattern = settings$time_point_pattern
)

# ===============================
# Build ROI Seurat objects
# ===============================

message("Building ROI Seurat objects...")

roi_objects <- build_roi_objects(
  roi_metadata = roi_metadata,
  seurat_obj = seurat_obj,
  sample_map = settings$sample_map
)

# ==============================
# Output directories
# ===============================

output_dirs <- make_output_dirs(
  analysis = "ROI_extraction",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "metadata",
    "objects"
  )
)

# ===============================
# Save ROI outputs
# ===============================

saveRDS(
  roi_metadata,
  file.path(
    output_dirs$metadata,
    paste0(
      settings$experiment_name,
    "_roi_metadata.rds"
    )
  )
)

saveRDS(
  roi_objects,
  file.path(
    output_dirs$objects,
    paste0(
      settings$experiment_name,
      "_roi_objects.rds"
    )
  )
)

# Done! 

message("Done!")
