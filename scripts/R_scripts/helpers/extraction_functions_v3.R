
# =====================================
# Function 1: Find ROI files
# =====================================

find_rois <- function(roi_dir) {
  
  fs::dir_ls(
    
    project_path(roi_dir),
    
    recurse = TRUE,
    regexp = "_cells_stats\\.csv$"
    
  )
  
}

# =====================================
# Function 2: Build ROI tibble
# =====================================

build_roi_tibble <- function(roi_files) {
  
  metadata <- settings$metadata
  
  roi_tibble <-
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
        settings$time_point_pattern
      ),
      
      # Folder above time_point
      experiment_folder = basename(dirname(dirname(file)))
      
    )
  
  # Check all time_points parsed correctly
  if (any(is.na(roi_tibble$time_point))) {
    
    stop(
      paste0(
        "Could not determine the time_point from:\n",
        paste(
          unique(dirname(roi_tibble$file[
            is.na(roi_tibble$time_point)
          ])),
          collapse = "\n"
        )
      )
      
    )
    
  }
  
  roi_tibble <-
    roi_tibble %>%
    separate(
      experiment_folder,
      into = metadata,
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
        group,
        time_point,
        sep = "_"
      )
      
    )
  
  roi_tibble <-
    roi_tibble %>%
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
  
  roi_tibble %>%
    relocate(
      experiment,
      group,
      tissue,
      time_point,
      roi,
      roi_id
    )
  
}

# =====================================
# Function 3: Create ROI Seurat object
# =====================================

create_roi_object <- function(
    roi_row,
    seurat_obj
) {
  
  roi_ids <- roi_row$cell_ids[[1]]
  
  matched_by_roi <- Cells(seurat_obj)[
    sub("_[0-9]+$", "", Cells(seurat_obj)) %in% roi_ids
  ]
  
  sample_counts <- table(seurat_obj$sample_id[matched_by_roi])
  
  sample_name <- names(which.max(sample_counts))
  
  matched_cells <- matched_by_roi[
    seurat_obj$sample_id[matched_by_roi] == sample_name
  ]
  
  if (max(sample_counts) / sum(sample_counts) < 0.95) {
    stop(
      "ROI is genuinely split across samples:\n",
      paste(names(sample_counts), sample_counts, collapse = ", ")
    )
  }
  
  if (length(matched_cells) == 0)
    stop(
      paste(
        "No matching cells found for",
        roi_row$roi_id
      )
    )
  
  if (length(unique(seurat_obj$sample_id[matched_cells])) > 1) {
    stop("Matched cells from multiple samples.")
  }
  
  cat(
    "ROI:",
    roi_row$group,
    roi_row$tissue,
    roi_row$time_point,
    roi_row$roi,
    "\n",
    
    "ROI cells:",
    length(roi_ids),
    "\n",
    
    "Matched:",
    length(matched_cells),
    "\n\n"
  )
  
  roi_obj <-
    subset(
      seurat_obj,
      cells = matched_cells
    )
  
  roi_obj$roi           <- roi_row$roi_id
  roi_obj$experiment    <- roi_row$experiment
  roi_obj$group         <- roi_row$group
  roi_obj$tissue        <- roi_row$tissue
  roi_obj$time_point     <- roi_row$time_point
  
  roi_obj
  
}

# =====================================
# Function 4: Build ROI objects
# =====================================

build_roi_objects <- function(
    roi_tibble,
    seurat_obj
) {
  
  roi_objects <-
    purrr::map(
      
      seq_len(nrow(roi_tibble)),
      
      ~ create_roi_object(
        roi_tibble[.x, ],
        seurat_obj
      )
      
    )
  
  names(roi_objects) <- roi_tibble$roi_id
  
  roi_objects

}

# =====================================
# Function 5: Extract ROIs
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