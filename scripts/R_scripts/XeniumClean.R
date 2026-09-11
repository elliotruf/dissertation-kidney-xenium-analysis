# ====================
# XeniumClean
# ====================
# For cleaning a Xenium dataset using XeniumClean.
# Xenium dataset must have its labels matched to a matching scRNA-seq dataset
# for the reference.
#
# - Set the label-transferred Xenium object for cleaning
#
#
#

# =====================
# User Settings
# =====================

settings <- list(
  
  xenium_obj = file.path(
    "results",
    "objects",
    "xen1_label_transfer.rds"
  ),
  
  reference_obj = file.path(
    "results",
    "objects",
    "xen1_sc_reference_new.rds"
  ),
  
  split_by = "time_point", # Most effective when the dataset is split by tissue slice, 
                          # i.e in xen1, by time point.
  
  new_object_name = "xen1_clean_bytimepoint_test2"
  
)

# =====================
# Dependencies
# =====================
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(tidyverse)
library(XeniumClean)

library(Banksy)

library(scater)
library(cowplot)
library(ggplot2)

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/celltype_functions_v3.R")

# ====================
# Load data
# ====================

xenium_lt <- readRDS(project_path(settings$xenium_obj))

reference <- readRDS(project_path(settings$reference_obj))

# =======================
# Prepare for cleaning
# =======================
# Ensure Assay5
reference[["RNA"]] <- as(reference[["RNA"]], Class = "Assay5")

message("Assigning broad cell types...")
# Create matching broad cell types
reference$broad_type <- assign_broad_celltypes(reference$cell_type)

xenium_lt$broad_type <- assign_broad_celltypes(xenium_lt$predicted.id)

message("Finding shared genes...")
# Find genes shared between datasets
shared_genes <- intersect(rownames(reference), rownames(xenium_lt))

# Split Xenium by "split_by" setting
xenium_list <- SplitObject(xenium_lt, split.by = settings$split_by)

# Join reference layers
reference_joined <- JoinLayers(
  reference,
  assay = "RNA"
)

message("Building reference gene set...")
# Build reference gene set
gene_sets <- XeniumClean::BuildReferenceGeneSets(
  reference = reference_joined,
  group.by = "broad_type",
  assay = "RNA",
  layer = "data",
  genes.use = shared_genes,
  expressed.threshold = 0.10,
  not.expressed.threshold = 0.05
)

message("Saving reference gene set...")
saveRDS(
  gene_sets,
  project_path(
    file.path(
      "results",
      "objects",
      "xen1_reference_gene_sets.rds"
    )
  )
)

# ===================
# Get FOV Coords
# ===================

xenium_list <- lapply(
  xenium_list,
  function(obj) {
    
    # Get all FOVs in the object
    fov_names <- names(obj@images)
    
    # Extract cell coordinates from each FOV
    coords_list <- lapply(
      fov_names,
      function(fov_name) {
        GetTissueCoordinates(
          obj,
          image = fov_name
        )
      }
    )
    
    # Combine coordinates from all FOVs
    coords <- do.call(rbind, coords_list)
    
    # Add coordinates to metadata using cell IDs
    obj$xenium_x <- coords$x[
      match(colnames(obj), coords$cell)
    ]
    
    obj$xenium_y <- coords$y[
      match(colnames(obj), coords$cell)
    ]
    
    obj
  }
)

# ===============================
# Clean dataset via XeniumClean
# ===============================
xenium_clean_list <- lapply(
  names(xenium_list),
  function(tp) {
    
    message("\n=== Cleaning ", tp, " ===")
    
    XeniumClean(
      object = xenium_list[[tp]],
      gene.sets = gene_sets,
      group.by = "broad_type",
      radius = 50,
      coords.source = "metadata",
      coord.cols = c("xenium_x", "xenium_y"),
      new.assay.name = "XeniumClean"
    )
  }
)

names(xenium_clean_list) <- names(xenium_list)

# ===========================
# Finalise cleaned dataset
# ===========================

message("Combining cleaned results...")
# Combine cleaned results
objs_bytimepoint <- unname(xenium_clean_list)

xenium_clean <- Reduce(
  function(x, y) merge(x, y, merge.data = TRUE),
  objs_bytimepoint
)

message("Running cleaned UMAP...")
# Run cleaned UMAP
DefaultAssay(xenium_clean) <- "XeniumClean"

# Normalise Cleaned Data
xenium_clean <- NormalizeData(xenium_clean, normalization.method = "RC", scale.factor = 100)
xenium_clean <- ScaleData(xenium_clean)

message("Saving cleaned Seurat object...")
# Save cleaned Seurat Object
saveRDS(
  xenium_clean,
  file = project_path(
    file.path(
      "results",
      "objects",
      paste0(settings$new_object_name, ".rds")
    )
  )
)
