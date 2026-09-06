# ================================================
# Functions for ROI Pseudobulk Analysis Pipeline
# ================================================

# ===============================
# Build Pseudobulk Profiles
# ===============================
# Sum gene counts across ROIs to generate pseudobulk
# count matrices and matching experiment metadata

build_pseudobulk <- function(
    roi_list,
    assay = "Xenium",
    layer = "counts"
) {
  
  counts <- lapply(
    roi_list,
    function(x) {
      
      Matrix::rowSums(
        SeuratObject::LayerData(
          x,
          assay = assay,
          layer = layer
        )
      )
      
    }
  )
  
  count_matrix <- do.call(cbind, counts)
  
  metadata <- purrr::map_dfr(
    roi_list,
    \(x) {
      
      tibble::tibble(
        roi = unique(x$roi),
        experiment = unique(x$experiment_name),
        condition = unique(x$condition),
        tissue = unique(x$tissue),
        time_point = unique(x$time_point)
      )
      
    }
  )
  
  metadata <- as.data.frame(metadata)
  
  metadata$sample_id <- make.unique(
    paste(
      metadata$roi,
      metadata$condition,
      metadata$time_point,
      sep = "_"
    )
  )
  
  rownames(metadata) <- metadata$sample_id
  colnames(count_matrix) <- metadata$sample_id
  
  message(
    "Count matrix dimensions: ",
    paste(dim(count_matrix), collapse = " x ")
  )
  
  message("Metadata rows: ", nrow(metadata))
  message("Sample IDs: ", length(metadata$sample_id))
  
  print(
    metadata[
      ,
      c("roi", "condition", "time_point", "sample_id")
    ]
  )
  
  metadata <- metadata[
    colnames(count_matrix),
    ,
    drop = FALSE
  ]
  
  list(
    counts = count_matrix,
    metadata = metadata
  )
}

# ==============================
# Prepare DE Metadata
# ==============================

prepare_de_metadata <- function(
    metadata,
    condition_variable,
    time_point_variable,
    reference_time_point
) {
  
  # Convert condition and time point to factors
  metadata[[condition_variable]] <-
    factor(metadata[[condition_variable]])
  
  metadata[[time_point_variable]] <-
    factor(metadata[[time_point_variable]])
  
  # Check that the reference time point exists
  if (!reference_time_point %in% levels(metadata[[time_point_variable]])) {
    stop(
      "Reference time point '",
      reference_time_point,
      "' not found in ",
      time_point_variable,
      "."
    )
  }
  
  # Set reference time point
  metadata[[time_point_variable]] <-
    relevel(
      metadata[[time_point_variable]],
      ref = reference_time_point
    )
  
  list(
    metadata = metadata
  )
}

# =======================================
# Fit edgeR Model
# =======================================
# Construct and fit an edgeR quasi-likelihood
# model from a pseudobulk count matrix.

fit_edge_r <- function(
    counts,
    metadata,
    variable
) {
  
  dge <- edgeR::DGEList(
    counts = counts,
    samples = metadata
  )
  
  # Filter lowly expressed genes.
  keep <- edgeR::filterByExpr(
    dge,
    group = metadata[[variable]]
  )
  
  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  if (nrow(dge) < 10) {
    stop("Too few genes remaining after filterByExpr.")
  }
  
  if (any(colSums(dge$counts) == 0)) {
    stop("One or more pseudobulks has zero counts after filtering.")
  }
  
  # Normalise library sizes.
  dge <- edgeR::calcNormFactors(dge)
  
  # Build design matrix.
  design <- model.matrix(
    as.formula(
      paste("~", variable)
    ),
    data = metadata
  )
  
  # Estimate dispersions.
  dge <- edgeR::estimateDisp(
    dge,
    design,
    robust = TRUE
  )
  
  # Fit quasi-likelihood model.
  fit <- edgeR::glmQLFit(
    dge,
    design,
    robust = TRUE
  )
  
  list(
    fit = fit,
    dge = dge,
    design = design
  )
}

# ========================================
# Run Differential Expression with EdgeR
# ========================================
# Run edgeR quasi-likelihood differential expression
# and return complete results table

run_roi_de <- function(
    fit,
    coef
) {
  
  qlf <- edgeR::glmQLFTest(
    fit,
    coef = coef
  )
  
  edgeR::topTags(
    qlf,
    n = Inf
  )$table
  
}

# ==========================================
# Run GO Enrichment
# ==========================================
# Perform GO enrichment for a single vector
# of gene symbols.

run_go <- function(
    genes,
    universe,
    ontology = "BP",
    p_cutoff = 0.05
) {
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  # Convert input genes to Entrez IDs.
  entrez <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # Convert Xenium panel genes to Entrez IDs.
  universe_entrez <- clusterProfiler::bitr(
    universe,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  if (is.null(entrez) || nrow(entrez) == 0) {
    return(NULL)
  }
  
  if (is.null(universe_entrez) || nrow(universe_entrez) == 0) {
    return(NULL)
  }
  
  clusterProfiler::enrichGO(
    gene = unique(entrez$ENTREZID),
    universe = unique(universe_entrez$ENTREZID),
    OrgDb = org.Mm.eg.db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff,
    readable = TRUE
  )
}

# ==========================================
# Run GO For All Comparisons
# ==========================================
# Perform GO enrichment separately for the
# up- and down-regulated genes from every
# differential expression comparison.

run_go_edge_r_all <- function(
    de_results,
    universe,
    fdr_cutoff,
    ontology,
    p_cutoff,
    min_genes
) {
  
  message(
    "Running GO analysis for ",
    length(de_results),
    " comparisons."
  )
  
  go_results <- list()
  
  for (comparison in names(de_results)) {
    
    sig <- subset(
      de_results[[comparison]],
      FDR < fdr_cutoff
    )
    
    up_genes <- rownames(
      subset(sig, logFC > 0)
    )
    
    down_genes <- rownames(
      subset(sig, logFC < 0)
    )
    
    # Upregulated genes.
    go_results[[paste0(comparison, "_up")]] <-
      if (length(up_genes) >= min_genes) {
        
        run_go(
          genes = up_genes,
          universe = universe,
          ontology = ontology,
          p_cutoff = p_cutoff
        )
        
      } else {
        
        NULL
        
      }
    
    # Downregulated genes.
    go_results[[paste0(comparison, "_down")]] <-
      if (length(down_genes) >= min_genes) {
        
        run_go(
          genes = down_genes,
          universe = universe,
          ontology = ontology,
          p_cutoff = p_cutoff
        )
        
      } else {
        
        NULL
        
      }
    
  }
  
  go_results
}
