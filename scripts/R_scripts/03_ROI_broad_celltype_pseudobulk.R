# Script for the analysis of ROIs at a broad_celltype level
# using pseudobulk + DE & GO

# ========================
# User Settings
# ========================

settings <- list(
  
  # ROI object to analyse
  roi_object = file.path(
    "results",
    "objects",
    "Xen1_Female_Cortex_roi_objects.rds"
  ),
  
  # Differential expression
  variable = "time_point",
  reference = "1wk",
  
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
# Preparation
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

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/pseudobulk_functions_v3.R")
source("scripts/R_scripts/helpers/broad_celltype_pseudobulk_functions_v2.R")
source("scripts/R_scripts/helpers/plotting_functions_v3.R")
source("scripts/R_scripts/helpers/output_functions_v3.R")

# ===============================
# Load ROI objects
# ===============================

message("Loading ROI object...")
rois <- readRDS(settings$roi_object)

run_name <- get_roi_run_name(rois)

output_dirs <- make_output_dirs(
  analysis = "ROI_broad_celltype_pseudobulk",
  run_name = run_name
)

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
# Build pseudobulks
# ===============================

message("Building pseudobulk...")
pb <- build_pseudobulk_broadtype(
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

message("Running differential expression analysis...")
de_results <- list()

# DE per broad celltype
for (broad in unique(metadata$broad_type)) {
  
  message("Running broad type: ", broad)
  
  # Metadata for this broad cell type
  metadata_ct <- subset(
    metadata,
    broad_type == broad
  )
  
  cat("\n========================\n")
  cat("Loop broad:", broad, "\n\n")
  
  print(unique(metadata_ct$broad_type))
  
  print(metadata_ct[, c("pseudobulk_id", "broad_type")])
  
  # Skip if too few pseudobulks
  if (nrow(metadata_ct) < 4) {
    message("Skipping ", broad, ": fewer than 4 ROI pseudobulks.")
    next
  }
  
  # Need two groups to compare
  if (length(unique(metadata_ct[[settings$variable]])) < 2)
    next
  
  # Counts for this broad cell type
  counts_ct <- counts[
    ,
    metadata_ct$pseudobulk_id,
    drop = FALSE
  ]
  
  # Need at least two replicates per group
  counts_per_group <- table(metadata_ct[[settings$variable]])
  
  if (any(counts_per_group < 2)) {
    message("Skipping", broad, ": fewer than two replicates.")
    next
  }
  
  fit_info <- tryCatch(
    fit_edgeR(
      counts = counts_ct,
      metadata = metadata_ct,
      variable = settings$variable
    ),
    error = function(e) {
      message("Skipping ", broad, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(fit_info))
    next
  
  reference <- levels(metadata_ct[[settings$variable]])[1]
  
  for (comp in comparisons) {
    
    coef_name <- paste0(
      settings$variable,
      comp
    )
    
    if (!coef_name %in% colnames(fit_info$design))
      next
    
    comparison_name <- paste(
      broad,
      paste0(comp, "_vs_", reference),
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
# Run GO enrichment
# ==========================================

message("Running Gene Ontology analysis...")

# Build Xenium panel gene universe
xenium_universe <- rownames(
  LayerData(
    rois[[1]],
    assay = settings$assay,
    layer = settings$layer
  )
)

# Run GO
go_results <- run_go_edger_all(
  de_results = de_results,
  universe = xenium_universe,
  fdr_cutoff = settings$fdr_cutoff,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)

message("Saving results...")
# =====================
# Save DE Tables
# =====================

save_de_tables_edgeR(
  de_results,
  output_dirs$de_tables,
  settings
)

# =====================
# Save DE Figures
# =====================

save_de_figures_edgeR(
  de_results,
  output_dirs$de_figures,
  settings,
  run_name
)

# =====================
# Save GO Tables
# =====================

save_go_tables(
  go_results,
  output_dirs$go_tables
)

# =====================
# Save GO Figures
# =====================

message("Saving GO figures")
save_go_figures_edgeR(
  go_results,
  output_dirs$go_figures,
  settings,
  run_name
)

# Done! 
message("Done!")