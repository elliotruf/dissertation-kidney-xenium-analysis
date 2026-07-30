
# Configurations for the processing and analysis of ROI selections

# Dataset Lookup
dataset_lookup <- list()

if (exists("xen1")) {
  dataset_lookup$xen1 <- list(
    seurat_obj = xen1,
    metadata = add_xen1_metadata
  )
}

if (exists("xen1diet")) {
  dataset_lookup$xen1diet <- list(
    seurat_obj = xen1diet,
    metadata = add_xen1_metadata
  )
}

if (exists("xen2")) {
  dataset_lookup$xen2 <- list(
    seurat_obj = xen2,
    metadata = add_xen2_metadata
  )
}

if (exists("xen2diet")) {
  dataset_lookup$xen2diet <- list(
    seurat_obj = xen2diet,
    metadata = add_xen2_metadata
  )
}
