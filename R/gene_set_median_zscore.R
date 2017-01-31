#' Determine median z-scores of a set of genes in a counts object
#' 
#' This function determines the patient-level median z-score of a set of genes. Counts of each gene are first converted to z-scores, then the median value for each patient is calculated.
#' @param gene_set the gene set, as a character vector.
#' @param counts a matrix or data frame of counts. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param remove_missing_genes logical, whether to remove genes from the gene sets if the genes are not present in the counts object. Defaults to TRUE.
#' @param remove_low_count_genes logical, whether to remove genes from the gene sets if the median expression level is below a threshold.  Implemented by a call to \code{gene_set_expressed}.
#' @param ... optional parameters passed to \code{gene_set_expressed}. Ignored if \code{remove_low_count_genes} is set to FALSE.
#' @export
#' @usage \code{gene_set_median_zscore(gene_set, counts,
#'   counts_genes_by="rows",
#'   remove_missing_genes=TRUE,
#'   remove_low_count_genes=TRUE,
#'   ...)}
#' @return A numeric vector containing the median values, with one element for each library in \code{counts}.
gene_set_median_zscore <-
  function(gene_set, counts,
           counts_genes_by="rows",
           remove_missing_genes=TRUE,
           remove_low_count_genes=TRUE,
           ...) {
    if (!is.character(gene_set))
      stop("Input object gene_set must be a character vector.")
    
    ## rotate matrix if needed
    counts_genes_by <- match.arg(counts_genes_by, choices=c("rows", "columns"))
    if (counts_genes_by=="columns") counts <- t(counts)
    
    ## remove missing genes if specified
    if (remove_missing_genes)
      gene_set <- gene_set[gene_set %in% rownames(counts)]
    
    if (remove_low_count_genes)
      gene_set <- gene_set_expressed(gene_set, counts, ...)
    
    ## stop if no genes remain
    if (length(gene_set) == 0)
      stop("No genes in gene set above expression threshold")
    
    ## convert counts to z-scores
    zscores <-
      scale(t(counts[match(gene_set, rownames(counts)), , drop=FALSE]))
    
    ## determine median value
    if (length(gene_set) == 1) {
      median_zscores <-
        zscores
    } else {
      median_zscores <-
        apply(zscores, MARGIN=1, median)
    }
    
    median_zscores
  }
