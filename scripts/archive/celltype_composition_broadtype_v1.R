# Script for the construction of a stacked barchart
# for the comparison of a single broad_type cell group

# =========================
# User Settings
# =========================

settings <- list(
  dataset_name = "xen1diet",
  region = "cortex",
  sex = "F",
  comparisons = c("1wk", "2wk", "4wk", "12wk"),
  # comparisons = c("Naive", "24h", "7d", "14d", "28d"),
  broad_type = "Stroma"
)

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

obj <- get(settings$dataset_name)

Idents(obj) <- "cell_type"

meta <- obj[[]]

meta$broad_type <- assign_broadtypes(meta$cell_type)

meta$time_point <- factor(
  meta$time_point,
  levels = settings$comparisons
)

plot_data <- subset(
  meta,
  broad_type == settings$broad_type &
  !is.na(time_point)
)

p_celltypeComposition <- ggplot(
  plot_data,
  aes(time_point, fill = cell_type)
) +
  
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = percent_format()
  ) +
  
  labs(
    x = "Developmental stage",
    y = "Cell Proportion",
    fill = "Cell subtype",
    title = "Cell Composition",
    subtitle = paste0(settings$dataset_name, "_", settings$broad_type)
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    "results",
    "composition",
    "figures",
    paste0(
      settings$dataset_name, "_",
      settings$region, "_",
      settings$broad_type,
      "_celltypeComposition.png"
    )
  ),
  
  plot = p_celltypeComposition,
  width = 6,
  height = 5
)

