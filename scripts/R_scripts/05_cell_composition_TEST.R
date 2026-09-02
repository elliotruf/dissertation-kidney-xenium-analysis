# Script for the construction of a stacked barchart
# for the comparison of a whole dataset and a broad cell type subset

# ========================
# User Settings
# ========================

settings <- list(
  
  # Experiment
  experiment = "xen1diet",
  
  # ROI object to analyse
  seurat_obj = file.path(
    "data",
    "xen1diet.rds"
  ),
  
  # Broad_celltype to analyse
  broad_celltype = "PT",
  
  # Condition column (sex, genotype, etc)
  condition = "sex", 
  
  # Order of time points for display
  comparisons = c("1wk", "2wk", "4wk", "12wk")
  # comparisons = c("Naive", "24h", "7d", "14d", "28d")
  
)

# =========================
# Preparation
# =========================

library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/celltype_functions_v3.R")
source("scripts/R_scripts/helpers/plotting_functions_v3.R")

# =============================
# Broad Cell Type Composition
# =============================

message("Loading data...")
seurat_obj <- readRDS(project_path(settings$seurat_obj))

message("Preparing object...")
Idents(seurat_obj) <- "cell_type"

seurat_obj <- add_broad_celltypes(seurat_obj)

meta <- seurat_obj[[]]

# Factor to order
meta$time_point <- factor(
  meta$time_point,
  levels = settings$comparisons
)

# Subset to specified broad_celltype
plot_data_broad <- subset(
  meta,
  broad_celltype == settings$broad_celltype &
    !is.na(time_point)
)

p_celltypeComposition_broad <- plot_cell_composition_broad(plot_data_broad)

# ====================================
# Whole Object Cell Type Composition
# ====================================

plot_data_whole <- subset(
  meta,
  !is.na(time_point)
) 

p_celltypeComposition_whole <- plot_cell_composition_whole(plot_data_whole)

# ===============
# Tables
# ===============

whole_summary <- summarise_whole_data(
  seurat_obj,
  condition_col = settings$condition,
  time_col = "time_point",
  broad_col = "broad_celltype",
  fine_col = "cell_type"
)

# ===============
# Save Outputs
# ===============

# Save tables
message("Saving tables...")
save_composition_summary_tables(
  whole_summary,
  project_path(
    "results",
    "cell_composition",
    "tables",
    settings$experiment
  )
)

# Save plots as pdf
message("Saving plots...")
ggsave(
  filename = file.path(
    "results",
    "cell_composition",
    "figures",
    paste0(
      settings$experiment, "_",
      settings$broad_celltype, "_",
      "celltype_composition.pdf"
    )
  ),
  
  plot = p_celltypeComposition_broad,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(
    "results",
    "cell_composition",
    "figures",
    paste0(
      settings$experiment,
      "_Whole_",
      "celltype_composition.pdf"
    )
  ),
  
  plot = p_celltypeComposition_whole,
  width = 6,
  height = 5
)

# Done!
message("Done!")