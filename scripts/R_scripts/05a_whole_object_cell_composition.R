# Script for the comparison of a Seurat object's cell type composition
# By fine- or broad-grained annotations.

# ========================
# User Settings
# ========================

settings <- list(
  
  # Dataset
  dataset_name = "xen2diet",
  
  # Seurat object
  seurat_obj = file.path(
    "data",
    "xen2diet.rds"
  ),
  
  # Experimental group
  group_column = "sample_id",
  
  # Cell type analysis level
  celltype_level = "broad",
  # "broad" = all broad cell types
  # "fine" = all fine cell types
  # "within_broad" = fine cell types within one broad type
  
  # Broad cell type to analyse when using "within_broad"
  parent_celltype = NULL,
  
  # Metadata columns
  fine_celltype_column = "cell_type",
  broad_celltype_column = "broad_celltype",
  timepoint_column = "time_point",
  
  # Order of time points for display
  timepoint_order = c(
    "sham",
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

source("scripts/R_scripts/helpers/project_paths_v3.R")
source("scripts/R_scripts/helpers/celltype_functions_v3.R")

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
# Assign Broad Cell Types
# ===============================

message("Assigning broad cell types...")

seurat_obj <- add_broad_celltypes(
  seurat_obj,
  fine_column = settings$fine_celltype_column
)

# ===============================
# Prepare data
# ===============================
meta <- seurat_obj[[]]

meta$time_point <- factor(
  meta[[settings$timepoint_column]],
  levels = settings$timepoint_order
)

meta <- subset(
  meta,
  !is.na(time_point)
)

# ===============================
# Check Required Columns
# ===============================

required_columns <- c(
  settings$group_column,
  settings$timepoint_column,
  settings$fine_celltype_column,
  settings$broad_celltype_column
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
  
}

# ===============================
# Select Cell Type Analysis
# ===============================

if (settings$celltype_level == "broad") {
  
  message("Analysing broad cell type composition...")
  
  meta$celltype_for_plot <-
    meta[[settings$broad_celltype_column]]
  
  plot_title <- "Broad cell-type composition"
  
  plot_subtitle <- settings$dataset_name
  
  legend_title <- "Broad cell type"
  
  output_name <- paste0(
    settings$dataset_name,
    "_broad_celltype_composition"
  )
  
}

if (settings$celltype_level == "fine") {
  
  message("Analysing fine cell type composition...")
  
  meta$celltype_for_plot <-
    meta[[settings$fine_celltype_column]]
  
  plot_title <- "Fine cell-type composition"
  
  plot_subtitle <- settings$dataset_name
  
  legend_title <- "Fine cell type"
  
  output_name <- paste0(
    settings$dataset_name,
    "_fine_celltype_composition"
  )
  
}

if (settings$celltype_level == "within_broad") {
  
  message(
    "Analysing fine cell type composition within ",
    settings$parent_celltype,
    "..."
  )
  
  meta <- subset(
    meta,
    meta[[settings$broad_celltype_column]] ==
      settings$parent_celltype
  )
  
  meta$celltype_for_plot <-
    meta[[settings$fine_celltype_column]]
  
  plot_title <- paste(
    "Fine cell-type composition within",
    settings$parent_celltype
  )
  
  plot_subtitle <- settings$dataset_name
  
  legend_title <- "Fine cell type"
  
  output_name <- paste0(
    settings$dataset_name,
    "_",
    paste(
      sort(unique(meta[[settings$group_column]])),
      collapse = "_"
    ),
    "_within_",
    settings$parent_celltype,
    "_fine_celltype_composition"
  )
  
}

# ===============================
# Calculate Cell Composition
# ===============================

plot_data <- meta %>%
  filter(
    !is.na(.data[[settings$group_column]]),
    !is.na(celltype_for_plot),
    celltype_for_plot != ""
  ) %>%
  count(
    .data[[settings$group_column]],
    time_point,
    celltype_for_plot,
    name = "cells"
  ) %>%
  group_by(
    .data[[settings$group_column]],
    time_point
  ) %>%
  mutate(
    proportion = cells / sum(cells)
  ) %>%
  ungroup()

# ===============================
# Plot Cell Composition
# ===============================

p_celltypeComposition <- ggplot(
  plot_data,
  aes(
    x = time_point,
    y = proportion,
    fill = celltype_for_plot
  )
) +
  
  geom_col() +
  
  facet_wrap(
    as.formula(
      paste("~", settings$group_column)
    )
  ) +
  
  scale_y_continuous(
    labels = percent_format()
  ) +
  
  labs(
    x = "Developmental stage",
    y = "Cell proportion",
    fill = legend_title,
    title = plot_title,
    subtitle = plot_subtitle
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

# ===============================
# Save Results
# ===============================

experiment_dir <- file.path(
  "results",
  "cell_composition",
  settings$dataset_name
)

figure_dir <- file.path(
  experiment_dir,
  "figures"
)

table_dir <- file.path(
  experiment_dir,
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
  width = 8,
  height = 9
)

# ===============================
# Save Composition Table
# ===============================

write.csv(
  plot_data,
  file = file.path(
    table_dir,
    paste0(
      output_name,
      ".csv"
    )
  ),
  row.names = FALSE
)

# Done!

message(
  "Done! Results saved using prefix: ",
  output_name
)