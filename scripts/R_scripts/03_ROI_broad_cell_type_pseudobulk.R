# ROI Broad Cell Type Pseudobulk
#
# Subsets each ROI into broad cell type groups,
# generates pseudobulk profiles, performs differential
# expression and Gene Ontology analyses, and saves
# volcano plots and tables.
#
# Before running:
# - Set the ROI object to the desired .rds file.
# - Set the experiment name used for output files
# - Set the condition, time point, and reference time point variables
# - Set the assay and layer
# - Set the desired FDR and LogFC thresholds, and the number of DEGs to be labeled.
# - Set the Gene Ontology to run (BP, MP, etc.), and the p-value threshold, 
# number of terms to display and the minimum required genes to run GO analysis.

# ========================
# User Settings
# ========================

settings <- list(
  
  # ROI object to analyse
  roi_object = file.path(
    "results",
    "ROI_extraction",
    "Xen1_Male_Cortex",
    "objects",
    "Xen1_Male_Cortex_roi_objects.rds"
  ),
  
  # Name used for output directories and files
  experiment_name = "Xen1_Male_Cortex",
  
  # Differential expression
  condition_variable = "condition",
  time_point_variable = "time_point",
  reference_time_point = "1wk",
  
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

# =========================
# Dependencies
# =========================

# Load modules
library(Seurat)
library(edgeR)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(pheatmap)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(enrichplot) 
library(stringr)

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/ROI_pseudobulk_functions.R")
source("scripts/R_scripts/helpers/ROI_broad_cell_type_pseudobulk_functions.R")
source("scripts/R_scripts/helpers/cell_type_functions.R")
source("scripts/R_scripts/helpers/plotting_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load ROI objects
# ===============================

message("Loading ROI object...")#

rois <- readRDS(settings$roi_object)

# ===============================
# Assign Broad Cell-Types
# ===============================

message("Assigning broad cell-types...")
rois <- lapply(
  rois,
  function(x) {
    
    x$broad_type <-
      assign_broad_celltypes(
        x$cell_type
      )
    
    x
    
  }
)

# ===============================
# Build Pseudobulk Profiles
# ===============================

message("Building pseudobulk profiles...")

pb <- build_pseudobulk_broadtype(
  roi_list = rois,
  assay = settings$assay,
  layer = settings$layer
)

counts <- pb$counts
metadata <- pb$metadata

# ===============================
# Prepare DE Metadata
# ===============================

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

message("Running differential expression analysis...")

de_results <- list()

# Run DE separately for each broad cell type.
for (broad_type in unique(metadata$broad_type)) {
  
  message("Running broad type: ", broad_type)
  
  # Metadata for this broad cell type.
  metadata_broad_type <- metadata[
    metadata$broad_type == broad_type,
    ,
    drop = FALSE
  ]
  
  cat("\n========================\n")
  cat("Broad type:", broad_type, "\n\n")
  
  print(
    metadata_broad_type[
      ,
      c("pseudobulk_id", "broad_type")
    ]
  )
  
  # Skip if too few pseudobulks.
  if (nrow(metadata_broad_type) < 4) {
    
    message(
      "Skipping ",
      broad_type,
      ": fewer than 4 ROI pseudobulks."
    )
    
    next
  }
  
  # Need at least two time points to compare.
  if (
    length(
      unique(
        metadata_broad_type[[settings$time_point_variable]]
      )
    ) < 2
  ) {
    
    message(
      "Skipping ",
      broad_type,
      ": fewer than two time points."
    )
    
    next
  }
  
  # Counts for this broad cell type.
  counts_broad_type <- counts[
    ,
    metadata_broad_type$pseudobulk_id,
    drop = FALSE
  ]
  
  # Require at least two pseudobulks per time point.
  counts_per_group <- table(
    metadata_broad_type[[settings$time_point_variable]]
  )
  
  if (any(counts_per_group < 2)) {
    
    message(
      "Skipping ",
      broad_type,
      ": fewer than two pseudobulks in one or more time points."
    )
    
    next
  }
  
  # Fit edgeR model.
  fit_info <- tryCatch(
    
    fit_edge_r(
      counts = counts_broad_type,
      metadata = metadata_broad_type,
      variable = settings$time_point_variable
    ),
    
    error = function(e) {
      
      message(
        "Skipping ",
        broad_type,
        ": ",
        e$message
      )
      
      NULL
      
    }
  )
  
  if (is.null(fit_info)) {
    next
  }
  
  reference_time_point <- levels(
    metadata_broad_type[[settings$time_point_variable]]
  )[1]
  
  comparisons <- setdiff(
    levels(
      metadata_broad_type[[settings$time_point_variable]]
    ),
    reference_time_point
  )
  
  for (comparison in comparisons) {
    
    coef_name <- paste0(
      settings$time_point_variable,
      comparison
    )
    
    if (!coef_name %in% colnames(fit_info$design)) {
      next
    }
    
    comparison_name <- paste(
      broad_type,
      paste0(
        comparison,
        "_vs_",
        reference_time_point
      ),
      sep = "_"
    )
    
    de_results[[comparison_name]] <-
      run_roi_de(
        fit = fit_info$fit,
        coef = coef_name
      )
    
  }
  
}

# ==========================================
# Run GO Enrichment
# ==========================================

message("Running Gene Ontology analysis...")

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

# ========================
# Make Output Directories
# ========================

output_dirs <- make_output_dirs(
  analysis = "ROI_broad_cell_type_pseudobulk",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "de_tables",
    "de_figures",
    "go_tables",
    "go_figures"
  )
)

# =========================
# Save Results
# =========================
message("Saving results...")

save_de_tables_edge_r(
  de_results = de_results,
  output_dir = output_dirs$de_tables,
  experiment_name = settings$experiment_name,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff
)

save_de_figures_edge_r(
  de_results = de_results,
  output_dir = output_dirs$de_figures,
  experiment_name = settings$experiment_name,
  fdr_cutoff = settings$fdr_cutoff,
  logfc_cutoff = settings$logfc_cutoff,
  n_labels = settings$n_labels
)

save_go_tables(
  go_results = go_results,
  output_dir = output_dirs$go_tables,
  experiment_name = settings$experiment_name
)

save_go_figures_edge_r(
  go_results = go_results,
  output_dir = output_dirs$go_figures,
  experiment_name = settings$experiment_name,
  n_terms = settings$go_n_terms
)

message("Done!")