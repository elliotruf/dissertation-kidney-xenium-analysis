# Functions for pseudobulk-based analyses

# Function 1: Pseudobulk ROIs (standard)
build_pseudobulk <- function(
    roi_list,
    assay = settings$assay,
    layer = settings$layer
) {
  
  # Build summed ROI matrix
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
  
  # Build metadata from ROI objects
  metadata <- purrr::map_dfr(
    roi_list,
    \(x) {
      data.frame(
        roi = dplyr::first(x$roi_id),
        sample = dplyr::first(x$sample),
        group = dplyr::first(x$group),
        tissue = dplyr::first(x$tissue),
        timepoint = dplyr::first(x$timepoint)
      )
    }
  )
  
  rownames(metadata) <- metadata$roi
  
  # Ensure metadata and counts are in same order
  metadata <- metadata[colnames(count_matrix), , drop = FALSE]
  
  list(
    counts = count_matrix,
    metadata = metadata
  )
}

# Function 2: Run DE with edgeR
run_de <- function(
    fit,
    coef,
    n = Inf
) {
  
  qlf <- glmQLFTest(
    fit,
    coef = coef
  )
  
  topTags(
    qlf,
    n = n,
  )$table
}

# Function 3: Run Gene Ontology analysis
run_go <- function(
    genes,
    ontology = settings$go_ontology
) {
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  # Convert symbol to ENTREZID
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # No mapped genes
  if (nrow(gene_df) == 0) {
    return(NULL)
  }
  
  enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = settings$go_p_cutoff,
    qvalueCutoff = settings$go_p_cutoff,
    readable = TRUE
  )
  
}

# Function 4: Naming convention
get_run_name <- function(rois) {
  
  sample <- unique(unlist(lapply(rois, \(x) x$sample)))
  group  <- unique(unlist(lapply(rois, \(x) x$group)))
  tissue <- unique(unlist(lapply(rois, \(x) x$tissue)))
  
  stopifnot(length(sample) == 1)
  stopifnot(length(group) == 1)
  stopifnot(length(tissue) == 1)
  
  paste(sample, group, tissue, sep = "_")
}
