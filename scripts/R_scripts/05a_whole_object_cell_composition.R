# ==========================================
# Whole-Kidney Cell Composition Analysis
# ==========================================
#
# Calculates cell type composition across time points for
# a selected experimental condition
# 
# Cell type analysis can be performed at:
#   - broad level (PT, Stroma, Vascular, Immune, etc.)
#   - fine level 
#   - fine level within one broad cell type (ex: all fine types within PT)
#
# Before running:
# - Set the experiment name.
# - Set the Seurat object.
# - Set the experimental conditions (column info is stored in, and desired condition)
# - Set the analysis level (broad, fine, within_broad)
# - Set the desired broad cell type (if using within_broad level)
# - Set the metadata columns storing broad and fine cell type columns, and time point
# - Set the time point order (Control, 24hrs, 48hrs, etc.)

# ========================
# User Settings
# ========================

settings <- list(
  
  # Name used for output directories and files
  experiment_name = "Xen2",
  
  # Seurat object
  seurat_obj = file.path(
    "data",
    "xen2diet.rds"
  ),
  
  # Experimental condition
  # Examples:
  #   Xen1: condition_variable = "sex"
  #         condition_value = "Female"
  #
  #   Xen2: condition_variable = "sample_id"
  #         condition_value = "wildType"
  condition_variable = "sample_id",
  condition_value = "wildType",
  
  # Cell-type analysis level
  #
  # "broad"       = all broad cell types
  # "fine"        = all fine cell types
  # "within_broad" = fine cell types within one broad type
  cell_type_level = "within_broad",
  
  # Broad cell type to analyse when using "within_broad"
  parent_cell_type = "Stroma",
  
  # Metadata columns
  fine_cell_type_column = "cell_type",
  broad_cell_type_column = "broad_cell_type",
  time_point_variable = "time_point",
  
  # Time-point order for display
  time_point_order = c(
    "Naive",
    "24h",
    "7d",
    "14d",
    "28d"
  )
)


# =========================
# Dependencies
# =========================

library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)

source("scripts/R_scripts/helpers/project_paths.R")
source("scripts/R_scripts/helpers/cell_type_functions.R")
source("scripts/R_scripts/helpers/output_functions.R")

# ===============================
# Load Seurat Object
# ===============================

message("Loading Seurat object...")

seurat_obj <- readRDS(
  project_path(settings$seurat_obj)
)

cat(
  "Loaded object with",
  ncol(seurat_obj),
  "cells.\n"
)

# ===============================
# Check required metadata
# ===============================

meta <- seurat_obj[[]]

required_columns <- c(
  settings$condition_variable,
  settings$time_point_variable,
  settings$fine_cell_type_column
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ========================
# Assign Broad Cell Types
# ========================

message("Assigning broad cell types...")

seurat_obj <- add_broad_cell_types(
  seurat_obj,
  fine_column = settings$fine_cell_type_column
)

meta <- seurat_obj[[]]

# ========================
# Check Condition
# ========================

available_conditions <- unique(
  meta[[settings$condition_variable]]
)

if (!settings$condition_value %in% available_conditions) {
  stop(
    "Condition value '",
    settings$condition_value,
    "' not found in ",
    settings$condition_variable,
    ". Available values are: ",
    paste(
      available_conditions,
      collapse = ", "
    )
  )
}

# ========================
# Prepare Seurat Object
# ========================

condition_cells <- seurat_obj[[]][[
  settings$condition_variable
]]

seurat_obj <- subset(
  seurat_obj,
  cells = rownames(seurat_obj[[]])[
    condition_cells == settings$condition_value
  ]
)

seurat_obj[[settings$time_point_variable]] <- factor(
  seurat_obj[[]][[
    settings$time_point_variable
  ]],
  levels = settings$time_point_order
)

# ========================
# Select Cell-Type Analysis
# ========================

if (settings$cell_type_level == "broad") {
  
  message(
    "Analysing broad cell-type composition..."
  )
  
  composition <- summarise_whole_data(
    seurat_obj = seurat_obj,
    condition_col = settings$condition_variable,
    time_col = settings$time_point_variable,
    broad_col = settings$broad_cell_type_column,
    fine_col = settings$fine_cell_type_column
  )
  
  plot_data <- composition$whole_composition_broad
  
  cell_type_column <- settings$broad_cell_type_column
  
  plot_data <- plot_data %>%
    filter(
      !is.na(.data[[cell_type_column]]),
      .data[[cell_type_column]] != ""
    )
  
  plot_title <- "Broad cell-type composition"
  
  legend_title <- "Broad cell type"
  
  output_name <- paste(
    settings$experiment_name,
    settings$condition_value,
    "broad_celltype_composition",
    sep = "_"
  )
  
} else if (settings$cell_type_level == "fine") {
  
  message(
    "Analysing fine cell-type composition..."
  )
  
  composition <- summarise_whole_data(
    seurat_obj = seurat_obj,
    condition_col = settings$condition_variable,
    time_col = settings$time_point_variable,
    broad_col = settings$broad_cell_type_column,
    fine_col = settings$fine_cell_type_column
  )
  
  plot_data <- composition$whole_composition_fine
  
  cell_type_column <- settings$fine_cell_type_column
  
  plot_data <- plot_data %>%
    filter(
      !is.na(.data[[cell_type_column]]),
      .data[[cell_type_column]] != ""
    )
  
  plot_title <- "Fine cell-type composition"
  
  legend_title <- "Fine cell type"
  
  output_name <- paste(
    settings$experiment_name,
    settings$condition_value,
    "fine_celltype_composition",
    sep = "_"
  )
  
} else if (settings$cell_type_level == "within_broad") {
  
  if (is.null(settings$parent_cell_type)) {
    stop(
      "parent_cell_type must be specified when ",
      "cell_type_level = 'within_broad'."
    )
  }
  
  message(
    "Analysing fine cell-type composition within ",
    settings$parent_cell_type,
    "..."
  )
  
  seurat_obj <- subset_cell_type(
    seurat_obj = seurat_obj,
    cell_type = settings$parent_cell_type,
    column = settings$broad_cell_type_column
  )
  
  composition <- summarise_whole_data(
    seurat_obj = seurat_obj,
    condition_col = settings$condition_variable,
    time_col = settings$time_point_variable,
    broad_col = settings$broad_cell_type_column,
    fine_col = settings$fine_cell_type_column
  )
  
  plot_data <- composition$whole_composition_fine
  
  cell_type_column <- settings$fine_cell_type_column
  
  plot_data <- plot_data %>%
    filter(
      !is.na(.data[[cell_type_column]]),
      .data[[cell_type_column]] != ""
    )
  
  plot_title <- paste(
    "Fine cell-type composition within",
    settings$parent_cell_type
  )
  
  legend_title <- "Fine cell type"
  
  output_name <- paste(
    settings$experiment_name,
    settings$condition_value,
    settings$parent_cell_type,
    "fine_celltype_composition",
    sep = "_"
  )
  
} else {
  
  stop(
    "Invalid cell_type_level: '",
    settings$cell_type_level,
    "'. Use 'broad', 'fine', or 'within_broad'."
  )
}


# ========================
# Plot Cell Composition
# ========================

p_celltype_composition <- ggplot(
  plot_data,
  aes(
    x = .data[[settings$time_point_variable]],
    y = proportion,
    fill = .data[[cell_type_column]]
  )
) +
  geom_col() +
  scale_x_discrete(
    limits = settings$time_point_order
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    x = "Time point",
    y = "Cell proportion",
    fill = legend_title,
    title = plot_title,
    subtitle = paste(
      settings$experiment_name,
      settings$condition_value,
      sep = " — "
    )
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  )

# ========================
# Create Output Directories
# ========================

output_dirs <- make_output_dirs(
  analysis = "cell_composition",
  experiment_name = settings$experiment_name,
  subdirectories = c(
    "figures",
    "tables"
  )
)

# ========================
# Save Plot
# ========================

ggsave(
  filename = file.path(
    output_dirs$figures,
    paste0(
      output_name,
      ".pdf"
    )
  ),
  plot = p_celltype_composition,
  width = 8,
  height = 9
)

# ========================
# Save Composition Table
# ========================

write.csv(
  plot_data,
  file = file.path(
    output_dirs$tables,
    paste0(
      output_name,
      ".csv"
    )
  ),
  row.names = FALSE
)

# Done!
message("Done!")