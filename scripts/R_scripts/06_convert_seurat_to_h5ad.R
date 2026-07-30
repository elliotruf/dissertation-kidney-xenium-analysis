# ==========================================
# Seurat -> .h5ad Object conversion
# ==========================================
# For converting a Seurat object to an .h5ad object
# in preparation for cell-cell communication inference using LIANA+ (Python).

# =========================================
# User Settings
# =========================================

settings <- list(
  
  experiment = "xen1diet",
  
  seurat_obj = file.path(
    "data",
    "xen1diet.rds"
  )
)

# =========================================

# Load dependencies
library(Seurat)
library(readxl)
library(SingleCellExperiment)
library(zellkonverter)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/celltype_functions_v3.R")
source("scripts/R_scripts/helpers/plotting_functions_v3.R")
source("scripts/R_scripts/helpers/pseudobulk_functions_v3.R")
source("scripts/R_scripts/helpers/wholekidney_celltype_functions_v2.R")

# Load data
message("Loading Seurat object...")
seurat_obj <- readRDS(project_path(settings$seurat_obj))

# Add broad celltype column and set assay
seurat_obj <- add_broad_celltypes(seurat_obj)

DefaultAssay(seurat_obj) <- "Xenium"

# Convert to .h5ad file
sce <- as.SingleCellExperiment(seurat_obj)

message("Saving .h5ad file...")
writeH5AD(
  sce,
  file = file.path("data", paste0(settings$experiment, ".h5ad"))
)

# Done!
message("Done!")
