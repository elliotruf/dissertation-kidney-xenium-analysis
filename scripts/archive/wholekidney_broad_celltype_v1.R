# ========================
# User Settings
# ========================

settings <- list(
  
  # Whole Seurat object to analyse
  seurat_obj = file.path(
    "xen1",
    "xen1diet.rds"
  ),
  
  broad_types = "Stroma",
  # Examples:
  # broad_types = "Stroma"
  # broad_types = c("Stroma", "Immune")
  # broad_types = NULL (to analyse ALL types)
  
  # Differential expression
  variable = "time_point",
  reference = "1wk",
  
  # edgeR
  assay = "Xenium",
  layer = "counts",
  
  # thresholds
  fdr_cutoff = 0.05,
  logfc_cutoff = 1,
  n_labels = 5,
  
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
library(dplyr)
library(Matrix)
library(pheatmap)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(enrichplot) 
library(stringr)

source("scripts/helpers/project_paths_v3.R")
source("scripts/helpers/wholekidney_broad_celltype_functions.R")
source("scripts/helpers/plotting_functions_v3.R")
source("scripts/helpers/output_functions_v3.R")

# ============================
# Load Seurat Object
# ============================

seurat_obj <- readRDS(settings$seurat_obj)

run_name <- get_run_name(
  seurat_obj,
  settings$seurat_obj
)

analysis_name <- if (is.null(settings$broad_types)) {
  
  "All"
  
} else {
  
  paste(settings$broad_types, collapse = "_")
  
}

output_dirs <- make_output_dirs(
  analysis = file.path(
    "wholekidney_broad_celltype",
    analysis_name
  ),
  run_name = run_name
)

# ===========================
# Assign Broad Cell Types
# ===========================

seurat_obj$broad_type <- assign_broad_celltypes(seurat_obj$cell_type)

# ==========================
# Select Broad Cell Types
# ==========================

broad_types <- if (is.null(settings$broad_types)) {
  
  sort(unique(seurat_obj$broad_type))
  
} else {
  
  settings$broad_types
  
}

missing_types <- setdiff(
  broad_types,
  unique(seurat_obj$broad_type)
)

if (length(missing_types) > 0) {
  stop(
    "Unknown broad type(s): ",
    paste(missing_types, collapse = ", ")
  )
}

message(
  "Running DE for: ",
  paste(broad_types, collapse = ", ")
)

# ==========================
# Run DE
# ==========================

de_results <- purrr::map(
  broad_types,
  ~run_broadtype_de(
    seurat_obj,
    .x,
    settings
  )
)

de_results <- purrr::list_flatten(de_results)

message("Finished DE")

# =====================
# Run GO
# =====================

go_results <- run_go_seurat(
  de_results = de_results,
  ontology = settings$go_ontology,
  p_cutoff = settings$go_p_cutoff,
  min_genes = settings$go_min_genes
)

# =====================
# Save DE Tables
# =====================

message("Saving DE tables")
save_de_tables_seurat(
  de_results,
  output_dirs$de_tables,
  settings
) 

# =====================
# Save DE Figures
# =====================

message("Saving volcanoes")
save_de_figures_seurat(
  de_results,
  output_dirs$de_figures,
  settings,
  run_name
)

# =====================
# Save GO Tables
# =====================

message("Saving GO tables")
save_go_tables(
  go_results,
  output_dirs$go_tables
)

# =====================
# Save GO Figures
# =====================

message("Saving GO figures")
save_go_figures_seurat(
  go_results,
  output_dirs$go_figures,
  settings,
  run_name
)

# Done! 
message("Done!")



