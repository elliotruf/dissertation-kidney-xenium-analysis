# DE and GO analysis according to broad cell type

# ===================================
# User Settings 
# ===================================

settings <- list(
  
  # Dataset
  experiment = "Xen1",
  
  # Seurat Object for analysis
  seurat_obj = file.path(
    "data",
    "xen1diet.rds"
  ),
  
  # Cell type for analysis
  celltype = "PT",
  # broad: "Stroma", "Immune", "PT", etc.
  # fine: 
  
  celltype_column = "broad_celltype",
  # "broad_celltype" for broad types
  # "cell_type" for fine types
  
  # Differential expression
  condition = "Female",
  variable = "time_point",
  reference = "1wk",
  time_point_order = c("1wk", "2wk", "4wk", "12wk"),
  
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

# ===============================
# Packages
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
source("scripts/R_scripts/helpers/wholekidney_celltype_functions.R")
source("scripts/R_scripts/helpers/plotting_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")
source("scripts/R_scripts/helpers/pseudobulk_functions.R")

# ===============================
# Prepare object
# ===============================

message("Loading data...")
seurat_obj <- readRDS(
  project_path(settings$seurat_obj)
)

cat("Class:", class(seurat_obj), "\n")

if (inherits(seurat_obj, "Seurat")) {
  
  cat(
    "Number of cells:",
    ncol(seurat_obj),
    "\n"
  )
  
} else {
  
  print(seurat_obj)
  
}

cat("Loaded object\n")

seurat_obj <- add_broad_celltypes(
  seurat_obj
)

cat("Added broad celltypes\n")

celltype_obj <- subset_celltype(
  seurat_obj,
  settings$celltype
)

celltype_obj <- subset(
  celltype_obj,
  subset = sample_id == settings$condition
)

cat(
  "Subset complete:",
  ncol(celltype_obj),
  "cells\n"
)

cat("Starting DE\n")

cell_summary <- summarise_celltype(
  celltype_obj,
  sample_col = "sample_id",
  variable = settings$variable
)

run_name <- get_celltype_run_name(
  settings
)

output_dirs <- make_output_dirs(
  analysis = "wholekidney_broad_celltype",
  run_name = run_name
)

# ===============================
# Differential Expression
# ===============================

message("Running differential expression analysis...")
celltype_obj[[settings$variable]] <- relevel(
  factor(celltype_obj[[settings$variable]][,1]), 
  ref = settings$reference
)

comparisons <- levels(celltype_obj[[settings$variable]][,1])

comparisons <- setdiff(
  comparisons,
  settings$reference
)

de_results <- list()

for (comp in comparisons) {
  
  comparison_name <- paste0(
    settings$condition,
    "_",
    comp,
    "_vs_",
    settings$reference
  )
  
  cat("Running:", comparison_name, "\n")
  
  de_results[[comparison_name]] <-
    run_celltype_de(
      celltype_obj,
      ident1 = comp,
      ident2 = settings$reference,
      variable = settings$variable,
      assay = settings$assay
    )
  
  cat("Finished:", comparison_name, "\n")
}

# ===============================
# Gene Ontology
# ===============================

message("Running GO analysis...")
# Build Xenium panel gene universe
xenium_universe <- rownames(
  LayerData(
    seurat_obj,
    assay = settings$assay,
    layer = settings$layer
  )
)

# Run GO (Seurat)
go_results <- run_all_go_seurat(
  de_results = de_results,
  universe = xenium_universe,
  fdr_cutoff = settings$fdr_cutoff,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)

# ===============================
# Figures
# ===============================

message("Creating figures...")
# Stacked DEGs bar plots
bar_plots <- lapply(
  names(de_results),
  function(name) {
    
    plot_top_genes_bar(
      de_table = de_results[[name]],
      title = name,
      run_name = run_name,
      fdr_cutoff = settings$fdr_cutoff,
      logfc_cutoff = settings$logfc_cutoff
    )
  }
)

# GO plots
go_plots <- purrr::imap(
  go_results,
  ~ plot_go_seurat(
    .x,
    title = .y
  )
)

# Plot DE Heatmap
heatmap <- plot_de_heatmap(
  de_results,
  seurat_obj = celltype_obj,
  group_by = settings$variable,
  assay = settings$assay,
  timepoint_order = settings$time_point_order
)

# =======================
# Save Outputs
# =======================

message("Saving results...")
# Save DE tables
save_de_tables_seurat(
  de_results,
  output_dirs$de_tables,
  settings
)

# Save DE figures
save_de_figures_seurat(
  de_results,
  output_dirs$de_figures,
  settings,
)

# Save GO tables
save_go_tables(
  go_results,
  output_dirs$go_tables
)

# Save GO figures
save_go_figures_seurat(
  go_results,
  output_dirs$go_figures,
  settings,
  run_name
)

# Save DE Heatmaps
save_de_heatmap_seurat(
  de_results,
  timepoint_order = settings$time_point_order,
  output_dir = output_dirs$de_figures,
  settings = settings,
  run_name = run_name,
  seurat_obj = celltype_obj,
  group_by = settings$variable,
  assay = settings$assay
)

# Done!
message("Done!")