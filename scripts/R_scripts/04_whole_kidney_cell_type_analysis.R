# ==========================================
# Whole Kidney Cell Type DE and GO Analysis
# ==========================================
#
# Subsets a Seurat object to a selected broad
# cell population and condition, performs
# differential expression across time points,
# runs Gene Ontology analysis, and saves
# tables and figures.
#
# Before running:
# - Set the Seurat object.
# - Set the cell population.
# - Set the condition.
# - Set the reference time point.
# - Set the time point order.
# - Set DE and GO thresholds.

# ===================================
# User Settings 
# ===================================

settings <- list(
  
  # Dataset
  experiment_name = "Xen1",
  
  # Seurat Object for analysis
  seurat_obj = file.path(
    "data",
    "xen1diet.rds"
  ),
  
  # Cell type for analysis
  # broad: "Stroma", "Immune", "PT", etc.
  # fine: 
  cell_type = "Immune",
  
  # Desired cell type variable
  #
  # "broad_cell_type" for broad cell types
  # "cell_type" for fine cell types
  cell_type_variable = "broad_cell_type",
 
  # Differential expression
  condition_variable = "sex",
  condition_value = "Female",
  time_point_variable = "time_point",
  reference_time_point = "1wk",
  time_point_order = c(
    "1wk",
    "2wk",
    "4wk",
    "12wk"
  ),
  
  # Assay
  assay = "Xenium",
  layer = "counts",
  
  # Thresholds
  fdr_cutoff = 0.05,
  logfc_cutoff = 1,
  n_labels = 10,
  
  # GO
  go_ontology = "BP",
  go_p_cutoff = 0.05,
  go_n_terms = 15,
  go_min_genes = 20
  
)

set.seed(1234)

# ===============================
# Dependencies
# ===============================

library(Seurat)
library(edgeR)
library(Matrix)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(presto)
library(readxl)

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/cell_type_functions.R")
source("scripts/R_scripts/helpers/ROI_pseudobulk_functions.R")
source("scripts/R_scripts/helpers/whole_kidney_cell_type_functions.R")
source("scripts/R_scripts/helpers/plotting_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load Data
# ===============================

message("Loading data...")

seurat_obj <- readRDS(
  project_path(settings$seurat_obj)
)

message(
  "Loaded object: ",
  class(seurat_obj)
)

if (inherits(seurat_obj, "Seurat")) {
  
  message(
    "Number of cells: ",
    ncol(seurat_obj)
  )
  
} else {
  
  print(seurat_obj)
  
}

# ===============================
# Assign Broad Cell Types
# ===============================

message("Assigning broad cell types...")

seurat_obj <- add_broad_cell_types(
  seurat_obj
)

message("Broad cell types assigned.")

# ===============================
# Subset Cell Population
# ===============================

message("Subsetting cell population...")

cell_type_obj <- subset_cell_type(
  seurat_obj,
  cell_type = settings$cell_type,
  column = settings$cell_type_variable
)

message(
  "Cells after cell population subset: ",
  ncol(cell_type_obj)
)

# ===============================
# Subset Condition
# ===============================

message("Subsetting condition...")

cell_condition <- cell_type_obj[[]][[settings$condition_variable]]

cell_type_obj <- subset(
  cell_type_obj,
  cells = colnames(cell_type_obj)[
    cell_condition == settings$condition_value
  ]
)

message(
  "Cells after condition subset: ",
  ncol(cell_type_obj)
)

# ===============================
# Differential Expression
# ===============================

message("Preparing time point metadata...")

cell_type_obj[[settings$time_point_variable]] <- relevel(
  factor(
    cell_type_obj[[settings$time_point_variable]][, 1],
    levels = settings$time_point_order
  ),
  ref = settings$reference_time_point
)

comparisons <- levels(
  cell_type_obj[[settings$time_point_variable]][, 1]
)

comparisons <- setdiff(
  comparisons,
  settings$reference_time_point
)

message(
  "Time point comparisons: ",
  paste(comparisons, collapse = ", ")
)

de_results <- list()

for (comp in comparisons) {
  
  comparison_name <- paste0(
    settings$cell_type,
    "_",
    settings$condition_value,
    "_",
    comp,
    "_vs_",
    settings$reference_time_point
  )
  
  message(
    "Running: ",
    comparison_name
  )
  
  de_results[[comparison_name]] <-
    run_cell_type_de(
      seurat_obj = cell_type_obj,
      group_1 = comp,
      group_2 = settings$reference_time_point,
      group_variable = settings$time_point_variable,
      assay = settings$assay
    )
  
  message(
    "Finished: ",
    comparison_name
  )
}

# ===============================
# Gene Ontology Analysis
# ===============================

message("Running GO analysis...")

# Build Xenium panel gene universe
xenium_universe <- rownames(
  SeuratObject::LayerData(
    seurat_obj,
    assay = settings$assay,
    layer = settings$layer
  )
)

go_results <- run_go_all_seurat(
  de_results = de_results,
  universe = xenium_universe,
  fdr_cutoff = settings$fdr_cutoff,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)

# ===============================
# Make Output Directories
# ===============================

output_dirs <- make_output_dirs(
  analysis = "whole_kidney_cell_type",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "de_tables",
    "de_figures",
    "go_tables",
    "go_figures"
  )
)

# ===============================
# Save Outputs
# ===============================

message("Saving results...")

# DE tables

save_de_tables_seurat(
  de_results = de_results,
  output_dir = output_dirs$de_tables,
  experiment_name = settings$experiment_name,
  condition_value = settings$condition_value,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff
)

# DE figures

save_de_figures_seurat(
  de_results = de_results,
  output_dir = output_dirs$de_figures,
  experiment_name = settings$experiment_name,
  condition_value = settings$condition_value,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff
)

# GO tables

save_go_tables(
  go_results = go_results,
  output_dir = output_dirs$go_tables,
  experiment_name = settings$experiment_name
)


# GO figures

save_go_figures_seurat(
  go_results = go_results,
  output_dir = output_dirs$go_figures,
  experiment_name = settings$experiment_name,
  condition_value = settings$condition_value,
  n_terms = settings$go_n_terms
)


# DE heatmap

save_de_heatmap_seurat(
  de_results = de_results,
  timepoint_order = settings$time_point_order,
  output_dir = output_dirs$de_figures,
  experiment_name = settings$experiment_name,
  cell_type = settings$cell_type,
  condition_value = settings$condition_value,
  seurat_obj = cell_type_obj,
  group_by = settings$time_point_variable,
  assay = settings$assay,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff
)

message("Done!")