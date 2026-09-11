# ======================================
# scRNA-seq to Xenium label transfer
# ======================================
# For the transfer of labels of an scRNAseq dataset to a Xenium dataset,
# in preparation for cleaning with XeniumClean (Zemek 2026).
#
# Before running:
# - Set the query Xenium object
# - Set the scRNA-seq reference object
# - Set the condition to split the reference by sample
# - Set the new object and reference names

# ===================
# User Settings
# ===================

settings <- list(
  
  # Xenium Seurat object
  xenium_obj = file.path(
    "xen1",
    "merged.samples.IndividualBanksy_AGFTrue_Lambda0.05_Resolution1.5_Annotated.rds"
  ),
  
  # scRNAseq Seurat object
  sc_obj = file.path(
    "data",
    "yuehan_stromal_integrated.rds"
  ),
  
  # Condition (to split scRNAseq set by sample)
  condition = "condition",
  
  new_object_name = "xen1_label_transfer",
  
  new_reference_name = "xen1_sc_reference_new"
  
)

# ========================
# Dependencies
# ========================

library(Seurat)
library(BPCells)
library(SeuratObject)
library(SeuratDisk)
library(tidyverse)
library(sctransform)
library(harmony)
library(jsonlite)
library(glmGamPoi)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(future)
library(conflicted)

source("scripts/R_scripts/helpers/project_paths_v3.R")

plan(sequential)

options(future.globals.maxSize = Inf)

# =======================
# Load Data
# =======================

message("Loading specified datasets...")

xenium_obj <- readRDS(project_path(settings$xenium_obj))

sc_obj <- readRDS(project_path(settings$sc_obj))

# ========================================
# Prepare sc_obj for label transfer
# ========================================
# Create new reference Seurat Object with gene symbols

# Extract raw counts
counts <- GetAssayData(
  sc_obj,
  assay = "RNA",
  layer = "counts"
)

# Convert Ensembl IDs to gene symbols
symbols <- mapIds(
  org.Mm.eg.db,
  keys = rownames(counts),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# Remove genes without symbols
keep <- !is.na(symbols)

counts <- counts[keep, ]
rownames(counts) <- symbols[keep]

# Remove duplicated gene symbols
counts <- counts[!duplicated(rownames(counts)), ]

# Create a new Seurat object
reference_new <- CreateSeuratObject(counts = counts)

# Copy metadata (cell annotations, clusters, etc.)
reference_new <- AddMetaData(
  reference_new,
  metadata = sc_obj@meta.data
)

# Normalize and compute PCA
reference_new <- NormalizeData(reference_new)
reference_new <- FindVariableFeatures(reference_new)
reference_new <- ScaleData(reference_new)
reference_new <- RunPCA(reference_new)

# Split reference by condition column 
sc.list <- SplitObject(
  reference_new,
  split.by = settings$condition
)

panel.genes <- sub(
  "\\.m[01]$",
  "",
  rownames(xenium_obj)
)

panel.genes <- unique(panel.genes)

residual_features <- intersect(
  panel.genes,
  rownames(reference_new)
)

reference.genes <- rownames(reference_new)

shared.genes <- intersect(
  panel.genes,
  reference.genes
)

# ============================================================================
# SCTransform and Harmony to prepare scRNA-seq for finding transfer anchors
# ============================================================================
# Group.by condition set to "condition"

# Transform each condition
sc.list <- lapply(
  X = sc.list,
  FUN = SCTransform,
  method = "glmGamPoi",
  return.only.var.genes = FALSE,
  residual.features = residual_features,
  verbose = FALSE
)

# Yield the union of all the features which were found to be variable 
# in at least one sc.list entry
var.features <- SelectIntegrationFeatures(
  object.list = sc.list,
  nfeatures = 3000
)

# List intersection of features
sc.sct <- merge(x = sc.list[[1]], y = sc.list[2:length(sc.list)], merge.data=TRUE) 

VariableFeatures(sc.sct) <- var.features

# Rescale only if residual features are supplied
if (!is.null(residual_features)){
  sc.sct <- ScaleData(sc.sct) # now, we will have all the genes in scale.data
}

sc.sct <- RunPCA(sc.sct, verbose = TRUE)

# Run Harmony by condition
sc.sct <- harmony::RunHarmony(
  object = sc.sct,
  group.by.vars = settings$condition
)

sc.sct <- RunUMAP(sc.sct, reduction = "harmony", dims = 1:30)

sc.sct <- FindNeighbors(sc.sct, reduction = "harmony", dims = 1:30) %>% FindClusters()

message("Single-cell object prepared.")

# ========================================================================
# SCTransform and Harmony to prepare Xenium for finding transfer anchors
# ========================================================================

message("Preparing Xenium object for label transfer...")
DefaultAssay(xenium_obj) <- "Xenium"

# Keep at least 6 transcripts in cell, preventing capture of bad quality or non-cells. 
xenium_obj <- subset(xenium_obj, subset = nCount_Xenium > 5)  
                                                            

xenium_obj <- SCTransform(
  xenium_obj, 
  assay = "Xenium", 
  vars.to.regress = c('nFeature_Xenium', 'nCount_Xenium'), 
  vst.flavor = "v2", 
  verbose = FALSE
)

xenium_obj <- RunPCA(xenium_obj, npcs = 50, features = rownames(data), verbose = FALSE)

xenium_obj <- RunUMAP(xenium_obj, dims = 1:30, verbose = FALSE)

xenium_obj <- FindNeighbors(xenium_obj, reduction = "pca", dims = 1:30, verbose = FALSE)

xenium_obj <- FindClusters(xenium_obj, resolution = 0.8, verbose = FALSE)

message("Xenium object prepared.")

# ========================
# Find Transfer Anchors
# ========================

message("Finding Transfer Anchors...")
anchors <- FindTransferAnchors(
  reference = sc.sct, 
  query = xenium_obj,
  normalization.method = "SCT",
  recompute.residuals = FALSE,
  reference.reduction = 'harmony'  # using batch corrected reduction
  ) 

message("Transferring labels...")
predictions <- TransferData(
  anchorset = anchors, 
  refdata = sc.sct$cell_type, 
  weight.reduction = xenium_obj[["pca"]], 
  dims = 1:20
  )

xenium_obj <- AddMetaData(
  xenium_obj,
  metadata = predictions
)

message("Saving new Xenium object with transferred labels...")
# Save new Xenium object with transferred label
saveRDS(
  xenium_obj,
  file = project_path(
    file.path(
      "results",
      "objects",
      paste0(settings$new_object_name, ".rds")
    )
  )
)

# Save reference with new SYMBOLs
saveRDS(
  sc.sct,
  file = project_path(
    file.path(
      "results",
      "objects",
      paste0(settings$new_reference_name, ".rds")
    )
  )
)

# Done!
message("Done!")

