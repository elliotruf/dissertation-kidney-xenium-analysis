# Script for the construction of a stacked barchart
# for the comparison of cell type

# =========================
# User Settings
# =========================

settings <- list(
  dataset_name = "xen1diet",
  region = "cortex",
  comparisons = c("1wk", "2wk", "4wk", "12wk"),
  # comparisons = c("Naive", "24h", "7d", "14d", "28d"),
)

library(Seurat)
library(ggplot2)
library(scales)

obj <- get(settings$dataset_name)

Idents(obj) <- "cell_type"

idents <- as.character(Idents(obj))

obj$broad_type <- assign_broadtype(ob)

obj$broad_type <- case_when(
  grepl("^PTS", idents) ~ "PT",
  grepl("^TAL", idents) ~ "TAL",
  grepl("^ATL|^DTL", idents) ~ "Thin limb",
  grepl("^CD", idents) ~ "Collecting duct",
  grepl("^DCT", idents) ~ "DCT",
  grepl("^STROMA", idents) ~ "Stroma",
  grepl("^ENDO", idents) ~ "Endothelium",
  grepl("^GEnC", idents) ~ "Endothelium",
  grepl("^VSMC", idents) ~ "Vascular",
  grepl("^IMMUNE", idents) ~ "Immune",
  grepl("^MESANGIAL", idents) ~ "Glomerular",
  grepl("^PODO", idents) ~ "Glomerular",
  grepl("^PEC", idents) ~ "Glomerular",
  grepl("^UROTHELIUM", idents) ~ "Urothelium",
  grepl("^FAT", idents) ~ "Adipose",
  grepl("^rare", idents) ~ "Rare",
  TRUE ~ "Other"
)

obj$time_point <- factor(
  obj$time_point,
  levels = settings$comparisons
)

plot_data <- subset(
  obj@meta.data,
  !is.na(time_point)
)

p_celltypeComposition <- ggplot(
  plot_data,
  aes(time_point, fill = broad_type$cell_type)
) +
  
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = percent_format()
  ) +
  
  labs(
    x = "Developmental stage",
    y = "Cell Proportion",
    fill = "Broad cell type",
    title = "Cell Type Composition",
    subtitle = settings$dataset
  )

ggsave(
  filename = file.path(
    "results",
    "composition",
    "figures",
    paste0(
      name, "_",
      settings$dataset_name, "_",
      settings$region,
      "_celltypeComposition.png"
    )
  ),
  
  plot = p_celltypeComposition,
  width = 6,
  height = 5,
)

