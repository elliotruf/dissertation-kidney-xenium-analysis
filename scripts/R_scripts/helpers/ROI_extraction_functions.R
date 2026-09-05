
# =====================================
# Find ROI files
# =====================================

find_roi_files <- function(roi_dir) {
  
  fs::dir_ls(
    
    project_path(roi_dir),
    
    recurse = TRUE,
    
    # Select only XeniumExplorer cell stats exports ("_cells_stats.csv")
    regexp = "_cells_stats\\.csv$"   
    
  )
  
}

# =====================================
# Build ROI metadata
# =====================================

build_roi_metadata <- function(
    roi_files,
    time_point_pattern
) {
  
  roi_metadata <-
    tibble(
      file = roi_files
    ) %>%
    mutate(
      
      filename = basename(file),
      
      # Folder containing ROI
      time_point = basename(dirname(file)),
      
      # Extract standardised time_point
      time_point = stringr::str_extract(
        time_point,
        time_point_pattern
      ),
      
      # Folder above time_point
      experiment_folder = basename(dirname(dirname(file)))
      
    )
  
  # Check all time_points parsed correctly
  if (any(is.na(roi_metadata$time_point))) {
    
    stop(
      paste0(
        "Could not determine the time_point from:\n",
        paste(
          unique(dirname(roi_metadata$file[
            is.na(roi$metadata_time_point)
          ])),
          collapse = "\n"
        )
      )
      
    )
    
  }
  
  roi_metadata <-
    roi_metadata %>%
    separate(
      experiment_folder,
      into = c("experiment_name", "condition", "tissue"),
      sep = "_",
      remove = FALSE,
      fill = "right"
    ) %>%
    mutate(
      
      roi = stringr::str_remove(
        filename,
        "_cells_stats\\.csv$"
      ),
      
      roi_id = paste(
        roi,
        condition,
        time_point,
        sep = "_"
      )
      
    )
  
  roi_metadata <-
    roi_metadata %>%
    mutate(
      
      roi_data = purrr::map(
        file,
        ~ readr::read_csv(
          .x,
          skip = 2,
          show_col_types = FALSE
        )
      ),
      
      cell_ids = purrr::map(
        roi_data,
        ~ .x[["Cell ID"]]
      )
      
    )
  
  roi_metadata %>%
    relocate(
      experiment_name,
      condition,
      tissue,
      time_point,
      roi,
      roi_id
    )
  
}

# =====================================
# Create ROI Seurat Object
# =====================================
# Create a single ROI object from one ROI row

create_roi_object <- function(
    roi_row,
    seurat_obj,
    sample_map = NULL
) {
  
  roi_ids <- roi_row$cell_ids[[1]]
  
  # Map ROI condition to Seurat sample ID
  sample_id <- if (is.null(sample_map)) {
    roi_row$condition
  } else {
    dplyr::recode(
      roi_row$condition,
      !!!sample_map
    )
  }
  
  # Restrict Seurat object to the expected sample first
  sample_cells <- Cells(seurat_obj)[
    seurat_obj$sample_id == sample_id
  ]
  
  if (length(sample_cells) == 0) {
    stop(
      "No Seurat cells found for sample: ",
      sample_id,
      " in ROI: ",
      roi_row$roi_id
    )
  }
  
  # Expected time point
  expected_timepoint <- tolower(roi_row$time_point)
  
  # Match ROI Cell IDs only within the expected sample
  matched_by_roi <- sample_cells[
    sub(
      "_[0-9]+$",
      "",
      sample_cells
    ) %in% roi_ids
  ]
  
  if (length(matched_by_roi) == 0) {
    cat("First few ROI IDs:\n")
    print(head(roi_ids))
    
    cat("First few sample cell names:\n")
    print(head(sample_cells))
    
    stop(
      "No matching cells found for ",
      roi_row$roi_id
    )
  }
  
  # Restrict matched cells to expected time point
  cell_timepoints <- tolower(
    seurat_obj$time_point[matched_by_roi]
  )
  
  cell_timepoints[cell_timepoints == "naive"] <- "sham"
  
  matched_by_roi <- matched_by_roi[
    cell_timepoints == expected_timepoint
  ]
  
  if (length(matched_by_roi) == 0) {
    stop(
      "No cells matching expected time point ",
      roi_row$time_point,
      " found for ",
      roi_row$roi_id
    )
  }
  
  # Strict sample QC
  matched_samples <- unique(
    seurat_obj$sample_id[matched_by_roi]
  )
  
  if (length(matched_samples) != 1 ||
      matched_samples != sample_id) {
    
    stop(
      "ROI matched unexpected sample(s): ",
      paste(matched_samples, collapse = ", ")
    )
  }
  
  # Strict time-point QC
  matched_timepoints <- unique(
    tolower(seurat_obj$time_point[matched_by_roi])
  )
  
  matched_timepoints[matched_timepoints == "naive"] <- "sham"
  
  if (length(matched_timepoints) != 1 ||
      matched_timepoints != expected_timepoint) {
    
    stop(
      "ROI ", roi_row$roi_id,
      " expected time point ", roi_row$time_point,
      " but matched cells from: ",
      paste(matched_timepoints, collapse = ", ")
    )
  }
  
  # Report matching
  # =====================================
  
  cat(
    "ROI:",
    roi_row$roi_id,
    "\n",
    
    "Expected sample:",
    sample_id,
    "\n",
    
    "Expected time point:",
    roi_row$time_point,
    "\n",
    
    "ROI cells:",
    length(roi_ids),
    "\n",
    
    "Matched:",
    length(matched_by_roi),
    "\n\n"
  )
  
  # ======================================
  
  # Create ROI object
  roi_obj <- subset(
    seurat_obj,
    cells = matched_by_roi
  )
  
  roi_obj$roi <- roi_row$roi_id
  roi_obj$experiment_name <- roi_row$experiment_name
  roi_obj$condition <- roi_row$condition
  roi_obj$tissue <- roi_row$tissue
  roi_obj$time_point <- roi_row$time_point
  
  roi_obj
  
}

# =====================================
# Build ROI objects
# =====================================
# Apply create_roi_object() to all ROI metadata
# and build the ROI collection

build_roi_objects <- function(
    roi_metadata,
    seurat_obj,
    sample_map = NULL
) {
  
  roi_objects <-
    purrr::map(
      
      seq_len(nrow(roi_metadata)),
      
      ~ create_roi_object(
        roi_row = roi_metadata[.x, ],
        seurat_obj = seurat_obj,
        sample_map = sample_map
      )
      
    )
  
  names(roi_objects) <- roi_metadata$roi_id
  
  roi_objects
  
}

# =====================================
# Extract ROIs
# =====================================

extract_rois <- function(seurat_obj) {
  
  roi_files <-
    find_rois(
      settings$roi_dir
    )
  
  roi_tibble <-
    build_roi_tibble(
      roi_files
    )
  
  roi_objects <-
    build_roi_objects(
      roi_tibble,
      seurat_obj
    )
  
  dir.create(
    project_path("results", "objects"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    
    roi_objects,
    
    file = project_path(
      
      "results",
      "objects",
      
      paste0(
        settings$experiment,
        "_roi_objects.rds"
      )
      
    )
    
  )
  
  saveRDS(
    
    roi_tibble,
    
    file = project_path(
      
      "results",
      "objects",
      
      paste0(
        settings$experiment,
        "_roi_tibble.rds"
      )
      
    )
    
  )
  
  invisible(
    
    list(
      
      roi_tibble = roi_tibble,
      roi_objects = roi_objects
      
    )
    
  )
  
}