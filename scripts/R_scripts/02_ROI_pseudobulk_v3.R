# ==================================
# User Settings
# ==================================

settings <- list(
  
  # ROI object to analyse
  roi_object = file.path(
    "results",
    "objects",
    "Xen2_WT_Vessels_roi_objects.rds"
  ),
  
  # Differential expression
  variable = "time_point",
  reference = "sham",
  
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
# Preparation
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

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/pseudobulk_functions_v3.R")
source("scripts/R_scripts/helpers/plotting_functions_v3.R")
source("scripts/R_scripts/helpers/output_functions_v3.R")

# ===============================
# Load ROI objects
# ===============================

rois <- readRDS(settings$roi_object)

run_name <- get_roi_run_name(rois)

output_dirs <- make_output_dirs(
  analysis = "ROI_pseudobulk",
  run_name = run_name
)

# ===============================
# Build pseudobulks
# ===============================

pb <- build_pseudobulk(
  rois,
  assay = settings$assay,
  layer = settings$layer
)

counts <- pb$counts
metadata <- pb$metadata

metadata_info <- prepare_de_metadata(
  metadata,
  variable = settings$variable,
  reference = settings$reference
)

metadata <- metadata_info$metadata

comparisons <- metadata_info$comparisons

# ==========================================
# Perform Differential Expression Analysis
# ==========================================

edgeR_model <- fit_edgeR(
  counts = counts,
  metadata = metadata,
  variable = settings$variable
)

fit <- edgeR_model$fit

reference <- levels(metadata[[settings$variable]])[1]

de_results <- list()

for (comp in comparisons) {
  
  coef_name <- paste0(
    settings$variable,
    comp
  )
  
  comparison_name <- paste0(
    comp,
    "_vs_",
    reference
  )
  
  de_results[[comparison_name]] <- run_roi_de(
    fit,
    coef_name
  )
  
}

# ===================================
# Perform Gene Ontology Analysis
# ===================================

go_results <- run_go_edger_all(
  de_results,
  fdr_cutoff = settings$fdr_cutoff,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)

print(names(go_results))

# ==========================================
# Save Differential Expression Results
# ==========================================

save_de_tables(
  de_results = de_results,
  output_dir = output_dirs$de_tables,
  settings = settings
)

# ==========================================
# Save Volcano Plots
# ==========================================

save_de_figures(
  de_results = de_results,
  output_dir = output_dirs$de_figures,
  settings = settings,
  run_name = run_name
)

# ==========================================
# Save Gene Ontology Results
# ==========================================

save_go_tables(
  go_results = go_results,
  output_dir = output_dirs$go_tables
)

# ==========================================
# Save Gene Ontology Figures
# ==========================================

save_go_figures(
  go_results = go_results,
  output_dir = output_dirs$go_figures,
  settings = settings,
  run_name = run_name
)