# ================================================
# Functions for ROI Pseudobulk Analysis Pipeline
# ================================================

# ===============================
# Function 1: Build Pseudobulks
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
        LayerData(
          x,
          assay = assay,
          layer = layer
        )
      )
      
    }
  )
  
  count_matrix <- do.call(cbind, counts)
  
  colnames(count_matrix) <- names(roi_list)
  
  metadata <- purrr::map_dfr(
    
    roi_list,
    
    \(x)
    
    tibble::tibble(
      
      roi = unique(x$roi),
      experiment = unique(x$experiment),
      group = unique(x$group),
      tissue = unique(x$tissue),
      time_point = unique(x$time_point)
      
    )
    
  )
  
  metadata <- as.data.frame(metadata)
  
  metadata$sample_id <- make.unique(
    paste(
      metadata$roi,
      metadata$group,
      metadata$time_point,
      sep = "_"
    )
  )
  
  # Use sample_id as common identifier
  rownames(metadata) <- metadata$sample_id
  
  # Give pseudobulk matrix the same column names
  colnames(count_matrix) <- metadata$sample_id
  
  # Reorder metadata to exactly match count matrix columns
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
# Function 2: Prepare DE Metadata
# ==============================

prepare_de_metadata <- function(
    metadata,
    variable,
    reference
) {
  
  metadata[[variable]] <-
    factor(metadata[[variable]])
  
  metadata[[variable]] <-
    relevel(
      metadata[[variable]],
      ref = reference
    )
  
  comparisons <-
    setdiff(
      levels(metadata[[variable]]),
      reference
    )
  
  list(
    metadata = metadata,
    comparisons = comparisons
  )
  
}


# =======================================
# Function 3: Fit edgeR Model
# =======================================
# Construct and fit an edgeR quasi-likelihood
# model from a pseudobulk count matrix

fit_edgeR <- function(
    counts,
    metadata,
    variable
) {
  
  # Build DGE object
  dge <- edgeR::DGEList(
    counts = counts,
    samples = metadata
  )
  
  # Filter lowly expressed genes
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
  
  # Normalise library sizes
  dge <- edgeR::calcNormFactors(dge)
  
  # Design matrix
  design <- model.matrix(
    as.formula(
      paste("~", variable)
    ),
    data = metadata
  )
  
  # Estimate dispersions
  dge <- edgeR::estimateDisp(
    dge,
    design,
    robust = TRUE
  )
  
  # Fit quasi-likelihood model
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
# Function 4: Run Differential Expression
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
# Function 3a: Run GO Enrichment
# ==========================================
# Perform GO enrichment for a single vector
# of gene symbols.

run_go <- function(
    genes,
    universe,
    ontology = "BP",
    p_cutoff = 0.05
) {
  
  if (length(genes) == 0)
    return(NULL)
  
  # Convert input genes to Entrez IDs
  entrez <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # Convert Xenium panel genes to Entrez IDs
  universe_entrez <- clusterProfiler::bitr(
    universe,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  if (is.null(entrez) || nrow(entrez) == 0)
    return(NULL)
  
  if (is.null(universe_entrez) || nrow(universe_entrez) == 0)
    return(NULL)
  
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
# Function 3b: Run GO For All Comparisons
# ==========================================
# Perform GO enrichment separately for the
# up- and down-regulated genes from every
# differential expression comparison.

run_go_edger_all <- function(
    de_results,
    universe,
    fdr_cutoff = settings$fdr_cutoff,
    ontology = settings$go_ontology,
    p_cutoff = settings$go_p_cutoff,
    min_genes = settings$go_min_genes
) {
  
  message("Number of comparisons: ", length(de_results))
  print(names(de_results))
  
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
    
    # Upregulated genes
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
    
    # Downregulated genes
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

# =============================================
# Function 4: Run Name
# =============================================
# Generate a unique identifier for the current
# analysis from the ROI metadata

get_roi_run_name <- function(rois) {
  
  experiment <- unique(vapply(rois, \(x) unique(x$experiment), character(1)))
  group  <- unique(vapply(rois, \(x) unique(x$group), character(1)))
  tissue <- unique(vapply(rois, \(x) unique(x$tissue), character(1)))
  
  stopifnot(length(experiment) == 1)
  stopifnot(length(group) == 1)
  stopifnot(length(tissue) == 1)
  
  paste(experiment, group, tissue, sep = "_")
}
