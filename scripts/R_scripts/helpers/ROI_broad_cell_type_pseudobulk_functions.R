# ======================================================
# Functions for ROI-level broad cell type pseudobulk
# ======================================================

# ===========================================
# Build Broad Cell Type Pseudobulk Profiles
# ===========================================
# Subset each ROI by broad cell type and sum
# gene counts to generate ROI-level pseudobulk
# profiles for each broad cell type.

build_pseudobulk_broadtype <- function(
    roi_list,
    assay = "Xenium",
    layer = "counts"
) {
  
  counts_list <- list()
  metadata_list <- list()
  
  for (roi_name in names(roi_list)) {
    
    roi <- roi_list[[roi_name]]
    
    for (broad_type_name in unique(roi$broad_type)) {
      
      roi_broad_type <- subset(
        roi,
        subset = broad_type == broad_type_name
      )
      
      # Skip empty broad types.
      if (ncol(roi_broad_type) == 0) {
        next
      }
      
      pseudobulk_id <- paste(
        roi_name,
        broad_type_name,
        sep = "_"
      )
      
      counts_list[[pseudobulk_id]] <-
        Matrix::rowSums(
          SeuratObject::LayerData(
            roi_broad_type,
            assay = assay,
            layer = layer
          )
        )
      
      # Check that each ROI/broad-type combination
      # contains cells from only one sample.
      if (length(unique(roi_broad_type$sample_id)) > 1) {
        
        message("Problem ROI: ", roi_name)
        message("Broad type: ", broad_type_name)
        
        print(
          unique(roi_broad_type$sample_id)
        )
        
        stop("Multiple sample_ids found.")
      }
      
      metadata_list[[pseudobulk_id]] <-
        tibble::tibble(
          pseudobulk_id = pseudobulk_id,
          experiment = unique(roi_broad_type$experiment_name),
          sample = unique(roi_broad_type$sample_id),
          roi = unique(roi_broad_type$roi),
          broad_type = broad_type_name,
          condition = unique(roi_broad_type$condition),
          tissue = unique(roi_broad_type$tissue),
          time_point = unique(roi_broad_type$time_point)
        )
      
    }
    
  }
  
  counts <- do.call(
    cbind,
    counts_list
  )
  
  metadata <- do.call(
    rbind,
    metadata_list
  )
  
  metadata <- as.data.frame(metadata)
  
  rownames(metadata) <- metadata$pseudobulk_id
  
  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  list(
    counts = counts,
    metadata = metadata
  )
}
