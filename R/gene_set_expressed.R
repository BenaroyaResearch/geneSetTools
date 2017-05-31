#' Trim low-count genes from gene sets
#' 
#' This function filters a gene set based on a set of criteria. The standard criterion is median gene expression above a certain threshold, but other functionality may be built in at a later time. It can optionally remove genes not present in the counts object, and/or remove genes from gene sets if their median counts are below a certain threshold.
#' @param gene_set the gene set, as a character vector.
#' @param counts a matrix or data frame of counts, or an object from which counts can be extracted. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param remove_missing_genes logical, whether to remove genes from the gene sets if the genes are not present in the counts object. Defaults to TRUE.
#' @param min_genes an integer, the minimum number of genes from a gene set that must be in the data. Gene sets with fewer genes overlapping with the dataset are dropped.
#' @param min_median_gene_set_expression numeric, the minimum value for median expression of genes in a gene set. Any gene sets with median expression lower than this value are removed.
#' @seealso \code{\link{filter_gene_sets}}
#' @import countSubsetNorm
#' @export
#' @usage \code{gene_set_expressed(gene_set, counts,
#'   counts_genes_by="rows",
#'   remove_missing_genes=TRUE,
#'   min_median_gene_expression=1)}
#' @return A character vector containing the genes that pass the criteria.  May be a zero-length vector, or identical to the input \code{gene_set} object.
gene_set_expressed <-
  function(gene_set, counts,
           counts_genes_by="rows",
           remove_missing_genes=TRUE,
           min_median_gene_expression=1) {
    if (!is.character(gene_set))
      stop("Input object gene_set must be a character vector.")
    
    counts <- extract_counts(counts)
    
    ## rotate matrix if needed
    counts_genes_by <- match.arg(counts_genes_by, choices=c("rows", "columns"))
    if (counts_genes_by=="columns") counts <- t(counts)
    
    ## remove missing genes if specified
    if (remove_missing_genes)
      gene_set <- gene_set[gene_set %in% rownames(counts)]
    
    ## determine median expression by gene
    if (length(gene_set) == 0) {
      return(gene_set)
    } else if (length(gene_set) == 1) {
      gene_set.median_exp <- median(counts[match(gene_set, rownames(counts)),])
    } else {
      gene_set.median_exp <-
        apply(counts[match(gene_set, rownames(counts)),],
              MARGIN=1, median)
    }
    
    ## remove genes not meeting expression threshold
    # keep genes not in counts object, if specified (by allowing NA median expression)
    gene_set <-
      na.omit(gene_set[is.na(gene_set.median_exp) |
                         gene_set.median_exp >= min_median_gene_expression])
    gene_set
  }
