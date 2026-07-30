# Functions for ROI extraction and analysis

# Function 1: Locate ROI files
find_rois <- function(roi_files) {
  
  roi_files <- list.files(
    path = settings$roi_dir,
    recursive = TRUE,
    pattern = "_cells_stats\\.csv$",
    full.names = TRUE
  )
}

# Function 2: Build ROI tibble
build_roi_tibble <- function(roi_files) {
  
  metadata <- settings$metadata
  
  roi_tibble <- tibble(
    file = roi_files,
  ) %>%
    mutate(
      filename = basename(file),
      
      ## Timepoints
      timepoint = basename(dirname(file)),
      
      ## Standardise Timepoints
      timepoint = str_extract(
        timepoint,
        settings$timepoint_pattern
      ),
      
      ## Sample
      sample_folder = basename(dirname(dirname(file)))
    )
  
  # Check all timepoints were parsed
  if (any(is.na(roi_tibble$timepoint))) {
    stop(
      "Could not determine the timepoint from the following folder(s):\n",
      paste(
        unique(
          basename(
            dirname(
              roi_tibble$file[is.na(roi_tibble$timepoint)]
        ))),
        collapse = ", "
      )
    )
  }
  
    # Split folder into parts
  roi_tibble <- roi_tibble %>%
    separate(
      sample_folder,
      into = metadata,
      sep = "_",
      remove = FALSE,
      fill = "right"
    ) %>%
    mutate(
      roi = str_remove(
        filename,
        "_cell_stats\\.csv$"
      ),
      roi_id = paste(
        sample,
        group,
        tissue,
        timepoint,
        roi,
        sep = "_"
      )
    )
  
  # Read ROI .csv
  roi_tibble <- roi_tibble %>%
    mutate(
      roi_data = purrr::map(
        file, 
        ~ read_csv(.x,
                   skip = 2,
                   show_col_types = FALSE)
      ),
      cell_ids = purrr::map(
        roi_data,
        ~ .x[["Cell ID"]]
      )
    )
  
  # Reorder columns by importance
  roi_tibble <- roi_tibble %>%
    relocate(
      sample,
      group,
      tissue,
      timepoint,
      roi
    )
  
  return(roi_tibble)
}

# Function 3: Create ROI Object
create_roi_object <- function(roi_row, seurat_obj) {
  
  roi_ids <- roi_row$cell_ids[[1]]
  
  matched_cells <- Cells(seurat_obj)[
    sub("_[0-9]+$", "", Cells(seurat_obj)) %in%
      roi_ids
  ]
  
  # Report matching
  if (length(matched_cells) == 0) {
    stop("No matching cells found.")
  }
  
  if (length(unique(matched_cells)) != length(matched_cells)) {
    warning("Duplicate matched cells detected")
  }
  
  cat(
    "ROI:", roi_row$group, roi_row$tissue, roi_row$timepoint, roi_row$roi, "\n",  
    "ROI cells:", length(roi_ids), "\n",
    "Matched:", length(matched_cells),"\n\n"
  )
  
  # Return ROI object
  roi_obj <- subset(seurat_obj, cells = matched_cells)
  
  roi_obj$roi <- roi_row$roi
  roi_obj$roi_id <- roi_row$roi_id
  roi_obj$sample <- roi_row$sample
  roi_obj$timepoint <- roi_row$timepoint
  roi_obj$tissue <- roi_row$tissue
  roi_obj$group <- roi_row$group
  
  return(roi_obj)
  
}

# Function 4: Build ROI Objects
build_roi_objects <- function(roi_tibble, seurat_obj) {
  
  roi_objects <- purrr::map(
    seq_len(nrow(roi_tibble)),
    ~ create_roi_object(
      roi_tibble[.x, ],
      seurat_obj
    )
  )
  
  names(roi_objects) <- roi_tibble$roi_id
  
  return(roi_objects)
}

# Function 5: Extract ROIs
extract_rois <- function(seurat_obj) {
  
  # Get ROIs from ROI directory
  roi_files <- find_rois(settings$roi_dir)
  
  # Build ROI tibble
  roi_tibble <- build_roi_tibble(roi_files)
  
  # Build ROI Seurat Objects
  roi_objects <- build_roi_objects(roi_tibble, seurat_obj)
  
  # Save ROI objects
  saveRDS(
    roi_objects,
    file = file.path(
      "results",
      "objects",
      paste0(settings$sample, "_roi_objects.rds")
    )
  )
  
  # Save ROI tibble
  saveRDS(
    roi_tibble,
    file = file.path(
      "results",
      "objects",
      paste0(settings$sample, "_roi_tibble.rds")
    )
  )
  
  # Return
  return(list(
    roi_tibble = roi_tibble,
    roi_objects = roi_objects
  ))
}
