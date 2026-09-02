# ======================================================
# Functions for ROI_broad_celltype_pseudobulk analysis
# ======================================================

# ==========================================
# Function 2: Build Broadtype Pseudobulks
# ==========================================

build_pseudobulk_broadtype <- function(
    roi_list,
    assay = "Xenium",
    layer = "counts"
) {
  
  counts_list <- list()
  metadata_list <- list()
  
  for (roi_name in names(roi_list)) {
    
    roi <- roi_list[[roi_name]]
    
    for (bt in unique(roi$broad_type)) {
      
      roi_bt <- subset(
        roi,
        subset = broad_type == bt
      )
      
      # Skip empty broadtypes
      if (ncol(roi_bt) == 0)
      next
      
      pseudobulk_id <- paste(
        roi_name,
        bt,
        sep = "_"
      )
      
      counts_list[[pseudobulk_id]] <-
        
        Matrix::rowSums(
          LayerData(
            roi_bt,
            assay = assay,
            layer = layer
          )
          
        )
      
      if (length(unique(roi_bt$sample_id)) > 1) {
        
        cat("\nProblem ROI:\n")
        print(roi_name)
        print(bt)
        
        print(unique(roi_bt$sample_id))
        
        stop("Multiple sample_ids found")
        
      }
      
      metadata_list[[pseudobulk_id]] <- 
        
        tibble::tibble(
          
          pseudobulk_id = pseudobulk_id,
          experiment = unique(roi_bt$experiment),
          sample = unique(roi_bt$sample_id),
          roi = unique(roi_bt$roi),
          broad_type = bt,
          group = unique(roi_bt$group),
          tissue = unique(roi_bt$tissue),
          time_point = unique(roi_bt$time_point),
          
          StringsAsFactors = FALSE
          
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
