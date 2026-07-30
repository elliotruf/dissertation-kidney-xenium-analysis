
# =========================
# User Settings
# =========================

settings <- list(
  roi_object = 
    "D:/USERS/ELLIOT/Dissertation/results/objects/Xen1_Male_Vessels_roi_objects.rds",
  
  variable = "timepoint",
  
  reference = "1wk", 
  
  # Analysis settings
  fdr_cutoff = 0.05, # Adjusted P-value (FDR) cut-off
  logfc_cutoff = 1, # LogFC cut-off
  n_labels = 10, # Top X DEGs
  
  assay = "Xenium",
  layer = "counts",

  go_ontology = "BP",
  go_p_cutoff = 0.05,
  go_n_terms = 15 # Top GO terms
)

# =========================
# Preparation
# =========================

# Load modules
library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(Matrix)
library(pheatmap)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(enrichplot) 
library(stringr)

source("scripts/helpers/pseudobulk_functions_v2.R")

# Extract ROIs
rois <- readRDS(settings$roi_object)

sample = dplyr::first(rois[[1]]$sample)
group = dplyr::first(rois[[1]]$group)
tissue = dplyr::first(rois[[1]]$tissue)

# =============================================
# Pseudobulk Differential Expression Analysis
# =============================================

# Build pseudobulks from ROIs
pb <- build_pseudobulk(rois)

counts <- pb$counts
metadata <- pb$metadata

# Ensure factor
metadata[[settings$variable]] <- factor(metadata[[settings$variable]])

# Set reference level
metadata[[settings$variable]] <- relevel(
  metadata[[settings$variable]],
  ref = settings$reference
)

reference <- levels(metadata[[settings$variable]])[1]

# Compare every other level against the reference
comparisons <- setdiff(
  levels(metadata[[settings$variable]]),
  reference
)

# Build DGE object with edgeR
dge <- DGEList(
  counts = counts,
  samples = metadata
)

# Filter genes by expression
keep <- filterByExpr(
  dge,
  group = metadata[[settings$variable]]
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

# Normalise
dge <- calcNormFactors(dge)

# Design matrix for edgeR
formula <- as.formula(
  paste("~", settings$variable)
)

design <- model.matrix(
  formula,
  data = metadata
)

# Estimate dispersal
dge <- estimateDisp(
  dge,
  design,
  robust = TRUE
)

# Fit model
fit <- glmQLFit(
  dge,
  design,
  robust = TRUE
)

# DE analysis
de_results <- list()

# Set naming convention
reference <- levels(metadata[[settings$variable]])[1]

# Run DE for each listed comparison
for (comp in comparisons) {
  
  coef_name <- paste0(
    settings$variable,
    comp
  )
  
  comparison_name <- paste0(
    comp,
    "_vs_",
    reference
  )
  
  de_results[[comparison_name]] <- run_de(
    fit,
    coef_name
  )
}

# ============================================
# Gene Ontology Analysis
# ============================================
go_results <- list()

for (name in names(de_results)) {
  
  sig_go <- subset(
    de_results[[name]],
    FDR < settings$fdr_cutoff
  )
  
  up_genes <- rownames(
    subset(sig_go, logFC > 0)
  )
  
  down_genes <- rownames(
    subset(sig_go, logFC < 0)
  )
  
  if (length(up_genes) >= 10) {
    go_results[[paste0(name, "_up")]] <- run_go(up_genes)
  } else {
    go_results[[paste0(name, "_up")]] <- NULL
  }
 
  if (length(down_genes) >= 10) {
    go_results[[paste0(name, "_down")]] <- run_go(down_genes)
  } else {
    go_results[[paste0(name, "_down")]] <- NULL
  }
  
}

# ============================================
# Save Results
# ============================================

# Naming convention
run_name <- get_run_name(rois)

output_dirs <- list(
  de_tables = file.path(
    "results",
    "ROI_pseudobulk",
    "DE",
    "tables",
    run_name
  ),
  
  go_tables = file.path(
    "results",
    "ROI_pseudobulk",
    "GO",
    "tables",
    run_name
  ),
  
  de_figures = file.path(
    "results",
    "ROI_pseudobulk",
    "DE",
    "figures",
    run_name
  ),
  
  go_figures = file.path(
    "results",
    "ROI_pseudobulk",
    "GO",
    "figures",
    run_name
  )
)

purrr::walk(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Save DE results
for (name in names(de_results)) {
  
  write.csv(
    de_results[[name]],
    file = file.path(
      output_dirs$de_tables,
      paste(
        name,
        "DE.csv",
        sep = "_")
      
      # Ex: "wk2_vs_wk1_xen1_cortex_F.csv" OR "WT_vs_KO_xen2_cortex.csv"
    )
  )
}

# Save significant genes only
for (name in names(de_results)) {
  
  sig_genes <- subset(
    de_results[[name]],
    FDR < settings$fdr_cutoff,
    abs(logFC) > settings$logfc_cutoff
  )

  write.csv(
    sig_genes,
    file = file.path(
      output_dirs$de_tables,
      paste(
        name,
        "DE_sig.csv",
        sep = "_")
      
      # Ex: "wk2_vs_wk1_xen1_cortex_F_sig.csv" OR "WT_vs_KO_xen2_cortex_sig.csv"
    )
  )
}

# Save volcano plots
for (name in names(de_results)) {
  
  # Get DE table
  df <- de_results[[name]]
  
  # Set significance
  df$Significance <- "Not significant"
  
  df$Significance[
    df$FDR < settings$fdr_cutoff & df$logFC > settings$logfc_cutoff
  ] <- "Up"
  
  df$Significance[
    df$FDR < settings$fdr_cutoff & df$logFC < -settings$logfc_cutoff
  ] <- "Down"

  # Significant genes only
  sig_genes <- subset(
    df,
    FDR < settings$fdr_cutoff &
      abs(logFC) > settings$logfc_cutoff
  )
  
  # Top genes to label
  if (nrow(sig_genes) > 0) {
    
  n_labels <- min(
    settings$n_labels,
    nrow(sig_genes)
  )
  
  top_sig <- sig_genes[
    order(sig_genes$FDR),
    ,
    drop = FALSE
  ]
  
  top_sig <- top_sig[
    seq_len(n_labels),
    ,
    drop = FALSE
  ]
  
  } else {
  
    top_sig <- NULL
}
  
  # Construct plot
  p_volcano <- ggplot(
    df,
    aes(
      x = logFC,
      y = -log10(FDR),
      color = Significance
    )
  ) +
  
    geom_point(alpha = 0.7) +
    
    scale_color_manual(
      values = c(
        "Up" = "red",
        "Down" = "blue",
        "Not significant" = "grey"
      )
    ) +
  
    geom_vline(
      xintercept = c(-1, 1),
      linetype = 2
    ) +
  
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2
    ) + 

    labs(
      title = name,
      x = "log2 Fold Change",
      y = "-log10 FDR"
    ) +
      theme_minimal()
  
  # Add labels to significant genes
  if (!is.null(top_sig)) {
    
    top_sig <- cbind(
      top_sig,
      Gene = rownames(top_sig)
    )
    
    p_volcano <- p_volcano +
      geom_text_repel(
        data = top_sig,
        aes(
          x = logFC,
          y = -log10(FDR),
          label = Gene
        ),
        
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = Inf
      )
  }
  
  ggsave(
    filename = file.path(
      output_dirs$de_figures,
      paste(
        name,
        "volcano.pdf",
        sep = "_"
      )
    ),
    
    plot = p_volcano,
    width = 6,
    height = 5,
  )
}

# Save GO analysis
for (name in names(go_results)) {
  
  write.csv(
    as.data.frame(go_results[[name]]),
    file.path(
      output_dirs$go_tables,
      paste(
        name,
        "GO.csv",
        sep = "_")
    ),
    row.names = FALSE
  )
  
  p_go <- enrichplot::dotplot(
    go_results[[name]],
    showCategory = settings$go_n_terms
  ) +
    scale_y_discrete(
      labels= \(x) str_wrap (x, width = 25)
    )
  
ggsave(
  filename = file.path(
    output_dirs$go_figures,
    paste(
      name,
      "GO.pdf",
      sep = "_"
    )
  ),
  
  plot = p_go,
  width = 10,
  height = 7
)

}