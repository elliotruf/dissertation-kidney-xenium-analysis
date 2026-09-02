# Script for the comparison of an ROI list's cell type composition
# By fine- or broad-grained annotations.

# ========================
# User Settings
# ========================

settings <- list(
  
  # Dataset
  dataset_name = "Xen2_KO_Vessels",
  
  # ROI object
  roi_object = file.path(
    "results",
    "objects",
    "Xen2_KO_Vessels_roi_objects.rds"
  ),
  
  # Cell type analysis level
  # "broad" = all broad cell types
  # "fine" = all fine cell types
  # "within_broad" = fine cell types within one broad type
  celltype_level = "broad",
  
  # Broad cell type to analyse when using "within_broad"
  parent_celltype = NULL,
  
  # Metadata columns
  fine_celltype_column = "cell_type",
  broad_celltype_column = "broad_type",
  timepoint_column = "time_point",
  
  # Order of time points
  timepoint_order = c(
    "Sham",
    "24h",
    "7d",
    "14d",
    "28d"
  )
  
)

# ===================
# Dependencies
# ===================

library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths_v3.R")

source("scripts/R_scripts/helpers/broad_celltype_pseudobulk_functions_v2.R")

# ===============================
# Load ROI Objects
# ===============================

message("Loading ROI objects...")

rois <- readRDS(
  project_path(settings$roi_object)
)

cat(
  "Loaded",
  length(rois),
  "ROIs.\n"
)

# ===============================
# Assign Broad Cell Types
# ===============================

# ===============================
# Assign Broad Cell Types
# ===============================

message("Assigning broad cell types...")

rois <- lapply(rois, function(roi) {
  
  roi$broad_type <- assign_broad_celltypes(
    roi[[settings$fine_celltype_column, drop = TRUE]]
  )
  
  roi$broad_type <- dplyr::case_when(
    roi$cell_type %in% c(
      "Endothelial_vessels"
    ) ~ "Endothelium",
    
    roi$cell_type %in% c(
      "PT_INJURY_acute_WT"
    ) ~ "PT",
    
    roi$cell_type %in% c(
      "CNT",
      "CNT _INJURY"
    ) ~ "DCT",
    
    roi$cell_type %in% c(
      "Podocytes_1",
      "Podocytes_2"
    ) ~ "Glomerular",
    
    roi$cell_type %in% c(
      "Urothelium"
    ) ~ "Urothelium",
    
    TRUE ~ roi$broad_type
  )
  
  roi
})

# ===============================
# Calculate ROI Composition
# ===============================

message("Calculating ROI cell composition...")

roi_composition_list <- lapply(
  names(rois),
  function(roi_name) {
    
    roi <- rois[[roi_name]]
    meta <- roi[[]]
    
    # Check time point column
    if (!settings$timepoint_column %in% colnames(meta)) {
      stop(
        "Time point column '",
        settings$timepoint_column,
        "' not found in ROI: ",
        roi_name
      )
    }
    
    # Get time point
    time_point <- unique(
      meta[[settings$timepoint_column]]
    )
    
    # Get expected time point from ROI name
    time_point <- sub(
      ".*_",
      "",
      roi_name
    )
    
    # Standardise Sham naming
    if (tolower(time_point) == "sham") {
      time_point <- "Sham"
    }
    
    if (is.na(time_point)) {
      stop(
        "Could not determine time point from ROI name: ",
        roi_name
      )
    }
    
    # Check the time points actually present
    actual_time_points <- unique(
      meta[[settings$timepoint_column]]
    )
    
    actual_time_points <- actual_time_points[
      !is.na(actual_time_points)
    ]
    
    # If multiple time points are present, retain only
    # cells matching the time point encoded in the ROI name
    if (length(actual_time_points) > 1) {
      
      message(
        "ROI ", roi_name,
        " contains multiple time points: ",
        paste(actual_time_points, collapse = ", "),
        ". Retaining only ",
        time_point,
        " cells."
      )
      
      keep <- meta[[settings$timepoint_column]] == time_point
      
      meta <- meta[
        keep,
        ,
        drop = FALSE
      ]
    }
    
    # Determine cell type column
    if (settings$celltype_level == "broad") {
      
      if (!settings$broad_celltype_column %in% colnames(meta)) {
        stop(
          "Broad cell type column '",
          settings$broad_celltype_column,
          "' not found in ROI: ",
          roi_name
        )
      }
      
      meta$celltype_for_plot <-
        meta[[settings$broad_celltype_column]]
      
    } else if (settings$celltype_level == "fine") {
      
      if (!settings$fine_celltype_column %in% colnames(meta)) {
        stop(
          "Fine cell type column '",
          settings$fine_celltype_column,
          "' not found in ROI: ",
          roi_name
        )
      }
      
      meta$celltype_for_plot <-
        meta[[settings$fine_celltype_column]]
      
    } else if (settings$celltype_level == "within_broad") {
      
      if (is.null(settings$parent_celltype)) {
        stop(
          "parent_celltype must be specified when ",
          "celltype_level = 'within_broad'."
        )
      }
      
      keep <- meta[[settings$broad_celltype_column]] ==
        settings$parent_celltype
      
      meta <- meta[
        keep,
        ,
        drop = FALSE
      ]
      
      meta$celltype_for_plot <-
        meta[[settings$fine_celltype_column]]
      
    } else {
      
      stop(
        "celltype_level must be 'broad', 'fine', ",
        "or 'within_broad'."
      )
      
    }
    
    # Remove missing annotations
    meta <- meta[
      !is.na(meta$celltype_for_plot) &
        meta$celltype_for_plot != "",
      ,
      drop = FALSE
    ]
    
    if (nrow(meta) == 0) {
      return(NULL)
    }
    
    # Count cells within each ROI
    result <- meta %>%
      count(
        celltype_for_plot,
        name = "cells"
      ) %>%
      mutate(
        roi = roi_name,
        time_point = time_point,
        proportion = cells / sum(cells)
      ) %>%
      select(
        roi,
        time_point,
        celltype = celltype_for_plot,
        cells,
        proportion
      )
    
    result
    
  }
)

roi_composition <- bind_rows(
  roi_composition_list
)

# ===============================
# Check Results
# ===============================

if (nrow(roi_composition) == 0) {
  stop(
    "No cell composition data were generated."
  )
}

message(
  "Calculated composition for ",
  length(unique(roi_composition$roi)),
  " ROIs."
)

# ===============================
# Prepare Time Point
# ===============================

roi_composition$time_point <- factor(
  roi_composition$time_point,
  levels = settings$timepoint_order
)

# ===============================
# Summarise Across ROIs
# ===============================

all_celltypes <- sort(unique(roi_composition$celltype))

plot_data <- roi_composition %>%
  group_by(time_point, roi) %>%
  tidyr::complete(
    celltype = all_celltypes,
    fill = list(proportion = 0)
  ) %>%
  ungroup() %>%
  group_by(time_point, celltype) %>%
  summarise(
    proportion = mean(proportion),
    .groups = "drop"
  )

# ===============================
# Plot Labels
# ===============================

if (settings$celltype_level == "broad") {
  
  legend_title <- "Broad cell type"
  
  plot_title <- paste(
    "Broad cell type composition:",
    settings$dataset_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
  
} else if (settings$celltype_level == "fine") {
  
  legend_title <- "Cell type"
  
  plot_title <- paste(
    "Cell type composition:",
    settings$dataset_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
  
} else {
  
  legend_title <- "Cell type"
  
  plot_title <- paste(
    "Cell composition within",
    settings$parent_celltype,
    ":",
    settings$dataset_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
  
}

# ===============================
# Plot Cell Composition
# ===============================

p_celltypeComposition <- ggplot(
  plot_data,
  aes(
    x = time_point,
    y = proportion,
    fill = celltype
  )
) +
  
  geom_col() +
  
  scale_y_continuous(
    labels = percent_format()
  ) +
  
  labs(
    x = "Time point",
    y = "Cell proportion",
    fill = legend_title,
    title = plot_title,
    subtitle = plot_subtitle
  ) +
  
  theme_classic()

print(
  p_celltypeComposition
)

# ===============================
# Output Name
# ===============================

if (settings$celltype_level == "broad") {
  
  output_name <- paste0(
    settings$dataset_name,
    "_ROI_broad_cell_composition"
  )
  
} else if (settings$celltype_level == "fine") {
  
  output_name <- paste0(
    settings$dataset_name,
    "_ROI_fine_cell_composition"
  )
  
} else {
  
  output_name <- paste0(
    settings$dataset_name,
    "_ROI_",
    settings$parent_celltype,
    "_fine_cell_composition"
  )
  
}

# ===============================
# Save Results
# ===============================

figure_dir <- file.path(
  "results",
  "cell_composition",
  "figures"
)

table_dir <- file.path(
  "results",
  "cell_composition",
  "tables"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ===============================
# Save Plot
# ===============================

ggsave(
  filename = file.path(
    figure_dir,
    paste0(
      output_name,
      ".pdf"
    )
  ),
  plot = p_celltypeComposition,
  width = 6,
  height = 5
)

# ===============================
# Save ROI-Level Composition
# ===============================

write.csv(
  roi_composition,
  file = file.path(
    table_dir,
    paste0(
      output_name,
      "_ROI_level.csv"
    )
  ),
  row.names = FALSE
)

# ===============================
# Save Summary
# ===============================

write.csv(
  plot_data,
  file = file.path(
    table_dir,
    paste0(
      output_name,
      "_summary.csv"
    )
  ),
  row.names = FALSE
)

# Done!
message(
  "Done! Results saved using prefix: ",
  output_name
)