# =============================================================
# Functions for the processing and analysis of ROI selections
# =============================================================

# Function 1: Extract & Subset ROIs from Seurat Object
extract_roi <- function(
    seurat_obj,
    roi_csv
) {
  
  # Read export
  roi <- read.csv(
    roi_csv,
    skip = 2
  )
  
  # Extract ROI IDs
  roi_ids <- roi$Cell.ID
  
  # Match ROI IDs to Seurat IDs
  cell_names <- rownames(seurat_obj@meta.data)
  
  matched_cells <- cell_names[
    sub("_[0-9]+$", "", cell_names) %in%
      roi_ids
  ]
  
  # Report matching
  cat(
    "ROI cells:",
    length(roi_ids),
    "\nMatched:",
    length(matched_cells),
    "\n"
  )
  
  if(length(matched_cells) == 0) {
    stop("No matching cells found.")
  }
  
  # Return ROI object
  subset(
    seurat_obj,
    cells = matched_cells
  )
}

# Function 2: Extract ROI files from folder and match
extract_roi_folder <- function(
    seurat_obj,
    folder_path
) {
  
  roi_files <- list.files(
    folder_path,
    pattern = "\\.csv$",
    full.names = TRUE
  )
  
  # Apply extract_roi() to each file in folder
  roi_list <- lapply(
    roi_files,
    function(x) {
      cat("Reading:", basename(x), "\n")
      extract_roi(seurat_obj, x)
  }
)
  
  # Remove "_cells_stats.csv" from file name
  roi_names <- gsub(
    "_cells_stats\\.csv$", 
    "",
    basename(roi_files)
  )
  
  names(roi_list) <- roi_names
  
  return(roi_list)
}

# Function 3: Build Cell Type Composition Dataframe
build_composition_df <- function(
    roi_list,
    celltype_column = "cell_type"
) {
  
  composition_df <- do.call(
    rbind,
    lapply(
      names(roi_list),
      function(roi) {
        
        comp <- prop.table(
          table(
            droplevels(
              roi_list[[roi]][[celltype_column]]
            )
          )
        ) * 100
        
        roi_df <- data.frame(
          ROI = roi,
          Cell_Type = names(comp), # Show selections in numerical order
          Percent = as.numeric(comp)
        )
        
        roi_df <- roi_df[
          order(-roi_df$Percent), # Show % by decreasing order
        ]
        
        return(roi_df)
      }
    )
  )
  
  rownames(composition_df) <- NULL
  
  return(composition_df)
}

# Function 4: Add Xen1 Metadata to Dataframe
add_xen1_metadata <- function(df) {
  
  df$ROI_Number <- as.numeric(
    sub(
      "Selection_([0-9]+).*",
      "\\1",
      df$ROI
    )
  )
  
  df$Sex <- sub(
    "Selection_[0-9]+_([FM])_.*",
    "\\1",
    df$ROI
  )
  
  df$Timepoint <- sub(
    ".*_(wk[0-9]+)$",
    "\\1",
    df$ROI
  )
  
  df$Timepoint <- factor(
    df$Timepoint,
    levels = c(
      "wk1",
      "wk2",
      "wk4",
      "wk12"
    )
  )
  
  print(class(df))
  str(df)
  
  df <- df[
    order(df$Timepoint, df$ROI_Number),
    ,
    drop = FALSE
  ]
  
  rownames(df) <- NULL
  
  return(df)
}

# Function 5: Add Xen2 Metadata to Dataframe

add_xen2_metadata <- NULL

# add_xen2_metadata <- function(df) {
  
#  df$ROI_Number <- as.numeric(
#    sub(
#      "Selection_([0-9]+).*",
#      "\\1",
#      df$ROI
#    )
#  )
    
#    df$Sex <- sub(
#      "Selection_[0-9]+_([FM])_.*",
#      "\\1",
#      df$ROI
#    )
    
#    df$Timepoint <- sub(
#      "Selection_"
#    )
    
#    df$Timepoint <- factor(
#      df$Timepoint,
#      levels = c(
#        "Naive",
#        "24h",
#        "7d",
#        "14d",
#        "28d"
#        ))
#}

# Function 6: Assign Broad Cell Types

assign_broadtypes <- function(cell_type) {
  
  case_when(
    grepl("^PTS", cell_type) ~ "PT",
    grepl("^TAL", cell_type) ~ "TAL",
    grepl("^ATL|^DTL", cell_type) ~ "Thin limb",
    grepl("^CD", cell_type) ~ "Collecting duct",
    grepl("^DCT", cell_type) ~ "DCT",
    grepl("^STROMA", cell_type) ~ "Stroma",
    grepl("^ENDO", cell_type) ~ "Endothelium",
    grepl("^GEnC", cell_type) ~ "Endothelium",
    grepl("^VSMC", cell_type) ~ "Vascular",
    grepl("^IMMUNE", cell_type) ~ "Immune",
    grepl("^MESANGIAL", cell_type) ~ "Glomerular",
    grepl("^PODO", cell_type) ~ "Glomerular",
    grepl("^PEC", cell_type) ~ "Glomerular",
    grepl("^UROTHELIUM", cell_type) ~ "Urothelium",
    grepl("^FAT", cell_type) ~ "Adipose",
    grepl("^rare", cell_type) ~ "Rare",
    TRUE ~ "Other"
  )
}

# Function 7: Pseudobulk ROIs (standard)
build_pseudobulk <- function(
    roi_list,
    metadata_fun = NULL,
    assay = "Xenium",
    layer = "counts"
) {
  
  # Build summed ROI matrix
  counts <- lapply(
    roi_list,
    function(x) {
      
      Matrix::rowSums(
        LayerData(
          x,
          assay = assay,
          layer = layer
        )
      )
    }
  )
  
  count_matrix <- do.call(
    cbind,
    counts,
  )
  
  colnames(count_matrix) <- names(roi_list)
  
  metadata <- NULL
  
  if (!is.null(metadata_fun)) {
    metadata <- data.frame(
      ROI = names(roi_list)
    )
    
    metadata <- metadata_fun(metadata)
    
    count_matrix <- count_matrix[
      ,
      match(metadata$ROI, colnames(count_matrix))
    ]
    
    return(list(
      counts = count_matrix,
      metadata = metadata
    ))
  }
}

# Function 8: Pseudobulk cell types across ROIs
build_pseudobulk_celltype <- function(
    roi_list,
    metadata_fun = NULL,
    assay = "Xenium",
    layer = "counts"
) {
  
  pseudobulks <- list()
  sample_info <- list()
  
  for (roi_name in names(roi_list)) {
    
    roi <- roi_list[[roi_name]]
    
    meta <- roi[[]]
  
    meta$broad_type <- assign_broadtypes(meta$cell_type)
    
    if (grepl("wk4", roi_name)) {
      cat("\n", roi_name, "\n")
      print(table(meta$broad_type))
    }
  
    for (broad in unique(meta$broad_type)) {
      
      cells <- rownames(
        meta[
          meta$broad_type == broad,
          ,
          drop = FALSE]
      )
      
      if (length(cells) == 0)
        next
        
      counts <- Matrix::rowSums(
        LayerData(
          roi,
          assay = assay,
          layer = layer
        )[, cells, drop = FALSE]
      )
      
      sample_name <- paste0(
        roi_name,
        "_",
        broad
      )
      
      pseudobulks[[sample_name]] <- counts
      
      sample_info[[sample_name]] <- data.frame(
        Sample = sample_name,
        ROI = roi_name,
        broad_type = broad
      )
    }
  }
   
    count_matrix <- do.call(
      cbind,
      pseudobulks
    )
    
    colnames(count_matrix) <- names(pseudobulks)
    
    metadata <- do.call(
      rbind,
      sample_info
    )
    
    # If metadata hasn't been previously added
    
    if (!is.null(metadata_fun)) {
      
      metadata <- metadata_fun(metadata)
    }
    
    return(
      list(
        counts = count_matrix,
        metadata = metadata
      )
    )
}

# Function 9: Run DE with edgeR
run_de <- function(
    fit,
    coef,
    n = Inf
) {
  
  qlf <- glmQLFTest(
    fit,
    coef = coef
  )
  
  topTags(
    qlf,
    n = n,
  )$table
}

# Function 10: Run Gene Ontology analysis
run_go <- function(
    genes,
    ontology = settings$go_ontology
) {
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  # Convert symbol to ENTREZID
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # No mapped genes
  if (nrow(gene_df) == 0) {
    return(NULL)
  }
  
  enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = settings$go_p_cutoff,
    qvalueCutoff = settings$go_p_cutoff,
    readable = TRUE
    )
  
}



