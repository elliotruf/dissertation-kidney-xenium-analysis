# ROI Pseudobulk
#
# Generates ROI-level pseudobulk profiles,
# performs differential expression and Gene Ontology analysis,
# and saves volcano plot visualisations and tables.
#
# Before running:
# - Set the ROI object to the desired .rds file.
# - Set the experiment name used for output files
# - Set the condition, time point, and reference time point variables
# - Set the assay and layer
# - Set the desired FDR and LogFC thresholds, and the number of DEGs to be labeled.
# - Set the Gene Ontology to run (BP, MP, etc.), and the p-value threshold, 
# number of terms to display and the minimum required genes to run GO analysis.

# ==================================
# User Settings
# ==================================

settings <- list(
  
  # ROI object to analyse
  roi_object = file.path(
    "results",
    "ROI_extraction",
    "Xen2_WT_Vessels",
    "objects",
    "Xen2_WT_Vessels_roi_objects.rds"
  ),
  
  # Name used for output directories and files
  experiment_name = "Xen2_WT_Vessels",
  
  # Differential expression
  condition_variable = "condition",   # What is the experimental condition variable called? 
  time_point_variable = "time_point", # What is the time point variable called?
  reference_time_point = "sham",      # What will the other time points be compared to?
  
  # edgeR
  assay = "Xenium",
  layer = "counts",
  
  # thresholds
  fdr_cutoff = 0.05,
  logfc_cutoff = 1,
  n_labels = 10,
  
  # GO
  go_ontology = "BP",
  go_p_cutoff = 0.05,
  go_n_terms = 15,
  go_min_genes = 5
  
)

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

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/ROI_pseudobulk_functions.R")
source("scripts/R_scripts/helpers/plotting_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load ROI objects
# ===============================

message("Loading ROI object...")
rois <- readRDS(settings$roi_object)

# Keep only cells matching the timepoint encoded in each ROI name
roi_names <- names(rois)

rois <- lapply(roi_names, function(roi_name) {
  
  x <- rois[[roi_name]]
  
  expected_timepoint <- sub(".*_", "", roi_name)
  
  subset(
    x,
    subset = time_point == expected_timepoint
  )
})

names(rois) <- roi_names

# ===============================
# Build Pseudobulk Profiles
# ===============================

message("Building pseudobulk...")

pb <- build_pseudobulk(
  roi_list = rois,
  assay = settings$assay,
  layer = settings$layer
)

counts <- pb$counts
metadata <- pb$metadata


# ==============================
# Prepare DE Metadata
# ==============================

metadata_info <- prepare_de_metadata(
  metadata = metadata,
  condition_variable = settings$condition_variable,
  time_point_variable = settings$time_point_variable,
  reference_time_point = settings$reference_time_point
)

metadata <- metadata_info$metadata


# ==========================================
# Perform Differential Expression Analysis
# ==========================================

message("Fitting edgeR model...")

edgeR_model <- fit_edge_r(
  counts = counts,
  metadata = metadata,
  variable = settings$time_point_variable
)

fit <- edgeR_model$fit

reference_time_point <- levels(
  metadata[[settings$time_point_variable]]
)[1]

comparisons <- setdiff(
  levels(metadata[[settings$time_point_variable]]),
  reference_time_point
)

de_results <- list()

message("Running differential expression analysis...")

for (comparison in comparisons) {
  
  coef_name <- paste0(
    settings$time_point_variable,
    comparison
  )
  
  comparison_name <- paste0(
    comparison,
    "_vs_",
    reference_time_point
  )
  
  de_results[[comparison_name]] <- run_roi_de(
    fit = fit,
    coef = coef_name
  )
  
}


# ===================================
# Perform Gene Ontology Analysis
# ===================================

message("Running Gene Ontology analysis...")

# Build Xenium panel gene universe.
xenium_universe <- rownames(
  SeuratObject::LayerData(
    rois[[1]],
    assay = settings$assay,
    layer = settings$layer
  )
)

go_results <- run_go_edge_r_all(
  de_results = de_results,
  universe = xenium_universe,
  fdr_cutoff = settings$fdr_cutoff,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)


# ===================================
# Make Output Directories
# ===================================

output_dirs <- make_output_dirs(
  analysis = "ROI_pseudobulk",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "de_tables",
    "de_figures",
    "go_tables",
    "go_figures"
  )
)

message("Saving results...")

# ==========================================
# Save DE Tables
# ==========================================

save_de_tables_edge_r(
  de_results = de_results,
  output_dir = output_dirs$de_tables,
  experiment_name = settings$experiment_name,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff
)

# ==========================================
# Save DE Figures
# ==========================================

save_de_figures_edge_r(
  de_results = de_results,
  output_dir = output_dirs$de_figures,
  experiment_name = settings$experiment_name,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff,
  n_labels = settings$n_labels
)


# ==========================================
# Save GO Tables
# ==========================================

save_go_tables(
  go_results = go_results,
  output_dir = output_dirs$go_tables,
  experiment_name = settings$experiment_name
)

# ==========================================
# Save GO Figures
# ==========================================

save_go_figures_edge_r(
  go_results = go_results,
  output_dir = output_dirs$go_figures,
  experiment_name = settings$experiment_name,
  n_terms = settings$go_n_terms
)


# Done
message("Done!")