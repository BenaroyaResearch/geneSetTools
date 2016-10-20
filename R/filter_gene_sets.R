#' Filter gene sets
#' 
#' This function filters one or more gene sets based on a set of criteria, removing gene sets that do not meet those criteria. The standard criterion is median gene expression above a certain threshold, but other functionality may be built in at a later time. It can optionally remove genes not present in the counts object, and/or remove genes from gene sets if their median counts are below a certain threshold.
#' @param gene_sets the gene set(s). Can be one of (1) a single gene set, as character vector, (2) a list of gene sets, with each as a character vector, or (3) a matrix or data frame, with one gene set per column, and the gene set names as column names. Names of list elements or columns are used as gene set names.
#' @param counts a matrix or data frame of counts. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param remove_missing_genes logical, whether to remove genes from the gene sets if the genes are not present in the counts object. Defaults to TRUE.
#' @param remove_low_count_genes logical, whether to remove genes from the gene sets if their median count does not exceed a threshold value. Defaults to FALSE.
#' @param min_median_gene_expression numeric, the minimum value for median expression of genes. Any genes with median expression lower than this value are removed from gene sets. Defaults to 1. Ignored if \code{remove_low_count_genes} is FALSE.
#' @param min_genes an integer, the minimum number of genes from a gene set that must be in the data. Gene sets with fewer genes overlapping with the dataset are dropped. Defaults to 2.
#' @param min_median_gene_set_expression numeric, the minimum value for median expression of genes in a gene set. Any gene sets with median expression lower than this value are removed.
#' @seealso \code{\link{gene_set_expressed}}
#' @export
#' @usage \code{filter_gene_sets(
#'              gene_sets, counts,
#'              counts_genes_by="rows",
#'              remove_missing_genes=TRUE,
#'              remove_low_count_genes=FALSE,
#'              min_median_gene_expression=1,
#'              min_genes=2,
#'              min_median_gene_set_expression=1)}
#' @return An object of the same class as the input object \code{gene_sets}, with the gene sets that meet the thresholds. 
filter_gene_sets <-
  function(gene_sets, counts,
           counts_genes_by="rows",
           remove_missing_genes=TRUE,
           remove_low_count_genes=FALSE,
           min_median_gene_expression=1,
           min_genes=2,
           min_median_gene_set_expression=1) {
    
    counts_genes_by <- match.arg(counts_genes_by, choices=c("rows", "columns"))
    if (counts_genes_by=="columns") counts <- t(counts)
    
    # store original class, and convert gene_sets to list (easier to do this and convert back at)
    class.orig <- class(gene_sets)
    if (is.data.frame(gene_sets)) {
      gene_sets <- as.list(gene_sets)
    } else if (is.matrix(gene_sets)) {
      gene_sets <- as.list(data.frame(gene_sets))
    } else if (!is.list(gene_sets)) {
      stop("gene_sets must be either a list, matrix, or data.frame.")
    }
    gene_sets <- lapply(gene_sets, na.omit) # remove NAs from each gene set
    
    # remove missing genes if requested
    if (remove_missing_genes)
      gene_sets <- lapply(gene_sets, function(x) x[x %in% rownames(counts)])
    
    ## remove low-count genes from gene sets, if requested
    # if there are 0-gene sets after this, they will be removed below
    if (remove_low_count_genes)
      gene_sets <- lapply(
        gene_sets, FUN=gene_set_expressed, counts=counts,
        remove_missing_genes=remove_missing_genes,
        min_median_gene_expression=min_median_gene_expression)
    
    ## calculate # of genes present
    gene_sets.ngenes_present <-
      sapply(gene_sets, function(x) sum(x %in% rownames(counts)))
    gene_sets.median_exp <-
      sapply(gene_sets, function(x) median(counts[na.omit(match(x, rownames(counts))),]))
    
    ## remove gene sets with too few genes or median expression too low
    gene_sets <-
      gene_sets[gene_sets.ngenes_present >= min_genes &
                  gene_sets.median_exp >= min_median_gene_set_expression]
    
    ## convert back to original class, and fill in matrix or data.frame with NAs
    if (class.orig %in% c("matrix", "data.frame")) {
      seq_max <- seq_len(max(sapply(gene_sets, length)))
      gene_sets <- sapply(gene_sets, "[", i=seq_max)
      if (class.orig=="data.frame")
        gene_sets <- data.frame(gene_sets)
    }
  
    gene_sets
  }
