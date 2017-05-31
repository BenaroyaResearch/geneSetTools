#' Determine median MAD-scaled values of a set of genes in a counts object
#' 
#' This function determines the patient-level median scaled value of a set of genes, centered on the median and scaled by the median absolute deviation (MAD). Counts of each gene are first converted to median-centered MAD-scaled values, then the median value for each patient is calculated.
#' @param gene_set the gene set, as a character vector.
#' @param counts a matrix or data frame of counts, or an object from which counts can be extracted. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param remove_missing_genes logical, whether to remove genes from the gene sets if the genes are not present in the counts object. Defaults to TRUE.
#' @param remove_low_count_genes logical, whether to remove genes from the gene sets if the median expression level is below a threshold.  Implemented by a call to \code{gene_set_expressed}.
#' @param ... optional parameters passed to \code{gene_set_expressed}. Ignored if \code{remove_low_count_genes} is set to FALSE.
#' @export
#' @usage \code{gene_set_median_mad_score(gene_set, counts,
#'   counts_genes_by="rows",
#'   remove_missing_genes=TRUE,
#'   remove_low_count_genes=TRUE,
#'   ...)}
#' @return A numeric vector containing the median values, with one element for each library in \code{counts}.
gene_set_median_mad_score <-
  function(gene_set, counts,
           counts_genes_by="rows",
           remove_missing_genes=TRUE,
           remove_low_count_genes=TRUE,
           ...) {
    if (!is.character(gene_set))
      stop("Input object gene_set must be a character vector.")
    
    counts <- extract_counts(counts)
    
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
    
    ## convert counts to mad scores
    mad_scores <-
      counts[match(gene_set, rownames(counts)), , drop=FALSE] %>%
      t() %>%
      apply(MARGIN=2, function(x) {(x - median(x))/mad(x)})
    
    ## determine median value
    if (length(gene_set) == 1) {
      median_mad_scores <-
        mad_scores
    } else {
      median_mad_scores <-
        apply(mad_scores, MARGIN=1, median)
    }
    
    median_mad_scores
  }
