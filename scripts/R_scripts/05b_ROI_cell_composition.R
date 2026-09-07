# ==========================================
# ROI-level Cell Composition Analysis
# ==========================================
#
# Calculates cell type composition across time points
# for a selected experimental condition.
#
# Cell type analysis can be performed at:
# - broad level
# - fine level
# - fine level within one broad cell type
#
# For "within_broad", only fine cell types belonging to
# the selected parent broad cell type are included.
#
# Before running:
# - Set the experiment name.
# - Set the ROI object.
# - Set the analysis level.
# - Set the desired broad cell type if using within_broad.
# - Set the metadata columns.
# - Set the time point order.

# ========================
# User Settings
# ========================

settings <- list(
  
  # Name used for output directories and files
  experiment_name = "Xen1_Male_Cortex",
  
  # ROI object
  roi_object = file.path(
    "results",
    "objects",
    "Xen1_Male_Cortex_roi_objects.rds"
  ),
  
  # Cell-type analysis level
  #
  # "broad"        = all broad cell types
  # "fine"         = all fine cell types
  # "within_broad" = fine cell types within one broad type
  cell_type_level = "broad",
  
  # Broad cell type to analyse when using "within_broad"
  parent_cell_type = NULL,
  
  # Metadata columns
  fine_cell_type_column = "cell_type",
  broad_cell_type_column = "broad_cell_type",
  time_point_variable = "time_point",
  
  # Time-point order for display
  time_point_order = c(
    "1wk",
    "2wk",
    "4wk",
    "12wk"
  )
)

# ===================
# Dependencies
# ===================

library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/cell_type_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load ROI Objects
# ===============================

message("Loading ROI objects...")
rois <- readRDS(
  project_path(settings$roi_object)
)

# ===============================
# Assign Broad Cell Types
# ===============================

message("Assigning broad cell types...")
rois <- lapply(
  rois,
  function(roi) {
    
    roi$broad_cell_type <- assign_broad_cell_types(
      roi[[settings$fine_cell_type_column, drop = TRUE]]
    )
    
    roi
    
  }
)

# ===============================
# Check Broad Cell Type Assignment
# ===============================

all_broad_cell_types <- unique(
  unlist(
    lapply(
      rois,
      function(roi) {
        roi[[settings$broad_cell_type_column, drop = TRUE]]
      }
    )
  )
)

# ===============================
# Validate Within-Broad Setting
# ===============================

if (settings$cell_type_level == "within_broad") {
  
  if (is.null(settings$parent_cell_type)) {
    stop(
      "parent_cell_type must be specified when ",
      "cell_type_level = 'within_broad'."
    )
  }
  
  if (!settings$parent_cell_type %in% all_broad_cell_types) {
    stop(
      "Selected parent cell type '",
      settings$parent_cell_type,
      "' was not found in broad cell type assignments."
    )
  }
}

# ===============================
# Calculate ROI Composition
# ===============================

message("Calculating ROI cell composition...")
roi_composition_list <- lapply(
  names(rois),
  function(roi_name) {
    
    roi <- rois[[roi_name]]
    meta <- roi[[]]
    
    # Check required columns
    
    required_columns <- c(
      settings$time_point_variable,
      settings$fine_cell_type_column,
      settings$broad_cell_type_column
    )
    
    missing_columns <- setdiff(
      required_columns,
      colnames(meta)
    )
    
    if (length(missing_columns) > 0) {
      stop(
        "Missing required column(s) in ROI ",
        roi_name,
        ": ",
        paste(missing_columns, collapse = ", ")
      )
    }
    
    # Get time point
    
    time_point <- unique(
      meta[[settings$time_point_variable]]
    )
    
    if (length(time_point) != 1) {
      stop(
        "Expected exactly one time point in ROI: ",
        roi_name,
        ". Found: ",
        paste(time_point, collapse = ", ")
      )
    }
    
    if (is.na(time_point)) {
      stop(
        "Time point is NA in ROI: ",
        roi_name
      )
    }
    
    # Select cell types to analyse
    
    if (settings$cell_type_level == "broad") {
      
      # Keep cells with a broad annotation
      meta <- meta[
        !is.na(meta[[settings$broad_cell_type_column]]),
        ,
        drop = FALSE
      ]
      
      # Plot broad cell types
      meta$cell_type_for_plot <- factor(
        meta[[settings$broad_cell_type_column]]
      )
      
    } else if (settings$cell_type_level == "fine") {
      
      # Keep cells with a fine annotation
      meta <- meta[
        !is.na(meta[[settings$fine_cell_type_column]]),
        ,
        drop = FALSE
      ]
      
      # Plot fine cell types
      meta$cell_type_for_plot <- factor(
        meta[[settings$fine_cell_type_column]]
      )
      
    } else if (settings$cell_type_level == "within_broad") {
      
      # Keep only the selected broad cell type
      keep <- !is.na(
        meta[[settings$broad_cell_type_column]]
      ) &
        meta[[settings$broad_cell_type_column]] ==
        settings$parent_cell_type
      
      meta <- meta[
        keep,
        ,
        drop = FALSE
      ]
      
      # Plot fine cell types within
      # the selected broad cell type
      meta$cell_type_for_plot <- factor(
        meta[[settings$fine_cell_type_column]]
      )
      
    } else {
      
      stop(
        "cell_type_level must be 'broad', 'fine', ",
        "or 'within_broad'."
      )
    }
    
    # Remove missing annotations
    
    meta <- meta[
      !is.na(meta$cell_type_for_plot) &
        meta$cell_type_for_plot != "",
      ,
      drop = FALSE
    ]
    
    # Drop unused factor levels after filtering
    meta$cell_type_for_plot <- droplevels(
      meta$cell_type_for_plot
    )
    
    if (nrow(meta) == 0) {
      warning(
        "No cells remaining for ",
        settings$cell_type_level,
        " analysis in ROI: ",
        roi_name
      )
      return(NULL)
    }
    
    # Count cells and calculate
    # proportions within this ROI
    
    meta %>%
      count(
        cell_type_for_plot,
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
        cell_type = cell_type_for_plot,
        cells,
        proportion
      )
    
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

# ===============================
# Prepare Time Point
# ===============================

roi_composition$time_point <- factor(
  roi_composition$time_point,
  levels = settings$time_point_order
)

# Stop if unexpected time points are present

if (any(is.na(roi_composition$time_point))) {
  stop(
    "One or more ROI time points were not found in ",
    "settings$time_point_order."
  )
}

# ===============================
# Prepare Cell-Type Factor
# ===============================

if (settings$cell_type_level == "broad") {
  
  roi_composition$cell_type <- factor(
    roi_composition$cell_type,
    levels = names(broad_cell_type_palette)
  )
  
  fill_palette <- broad_cell_type_palette
  
} else {
  
  roi_composition$cell_type <- factor(
    roi_composition$cell_type,
    levels = names(fine_cell_type_palette)
  )
  
  fill_palette <- fine_cell_type_palette
  
}
# ===============================
# Summarise Across ROIs
# ===============================

message("Summarising composition across ROIs...")
plot_data <- roi_composition %>%
  group_by(
    roi,
    time_point
  ) %>%
  tidyr::complete(
    cell_type,
    fill = list(
      cells = 0,
      proportion = 0
    )
  ) %>%
  ungroup() %>%
  group_by(
    time_point,
    cell_type
  ) %>%
  summarise(
    proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  filter(
    proportion > 0
  )

# ===============================
# Plot Labels
# ===============================

if (settings$cell_type_level == "broad") {
  
  legend_title <- "Broad cell type"
  
  plot_title <- paste(
    "Broad cell type composition:",
    settings$experiment_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
  
} else if (settings$cell_type_level == "fine") {
  
  legend_title <- "Cell type"
  
  plot_title <- paste(
    "Fine cell type composition:",
    settings$experiment_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
  
} else {
  
  legend_title <- "Fine cell type"
  
  plot_title <- paste(
    "Fine cell type composition within",
    settings$parent_cell_type,
    ":",
    settings$experiment_name
  )
  
  plot_subtitle <- "Mean composition across ROIs"
}

# ==============================
# Plot Cell Composition
# ==============================

p_cell_type_composition <- ggplot(
  plot_data,
  aes(
    x = .data[[settings$time_point_variable]],
    y = proportion,
    fill = cell_type
  )
) +
  geom_col() +
  scale_x_discrete(
    limits = settings$time_point_order
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = fill_palette,
    breaks = unique(plot_data$cell_type),
    drop = TRUE
  ) +
  labs(
    x = "Time point",
    y = "Cell proportion",
    fill = legend_title,
    title = plot_title,
    subtitle = plot_subtitle
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.spacing.y = unit(2, "pt"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm")
  ) +
  guides(
    fill = guide_legend(
      ncol = if (settings$cell_type_level == "broad") 3 else 4,
      byrow = TRUE
    )
  )

# ===============================
# Output Name
# ===============================

if (settings$cell_type_level == "broad") {
  
  analysis_name <- "broad_cell_type_composition"
  
} else if (settings$cell_type_level == "fine") {
  
  analysis_name <- "fine_cell_type_composition"
  
} else {
  
  analysis_name <- paste(
    "within_broad",
    settings$parent_cell_type,
    "cell_type_composition",
    sep = "_"
  )
}

output_name <- paste(
  settings$experiment_name,
  analysis_name,
  sep = "_"
)

# ===============================
# Make Output Directories
# ===============================

output_dirs <- make_output_dirs(
  analysis = "cell_composition",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "figures",
    "tables"
  )
)

# ===============================
# Save Plot
# ===============================

message("Saving results...")
ggsave(
  filename = file.path(
    output_dirs$figures,
    paste0(
      output_name,
      ".pdf"
    )
  ),
  plot = p_cell_type_composition,
  width = 7,
  height = 5
)

# ===============================
# Save ROI-Level Composition
# ===============================

readr::write_csv(
  roi_composition,
  file.path(
    output_dirs$tables,
    paste0(
      output_name,
      "_ROI_level.csv"
    )
  )
)

# ===============================
# Save Summary
# ===============================

readr::write_csv(
  plot_data,
  file.path(
    output_dirs$tables,
    paste0(
      output_name,
      "_summary.csv"
    )
  )
)
