#' Filter gene sets
#' 
#' This function filters a gene set based on a set of criteria. The standard criterion is median gene expression above a certain threshold, but other functionality may be built in at a later time. 
#' @param gene_sets the gene sets. Either a list of gene sets, with each as a character vector, or a matrix or data frame, with one gene set per column, and the gene set names as column names. Names of list elements are used as gene set names.
#' @param counts a matrix or data frame of counts. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param min_genes an integer, the minimum number of genes from a gene set that must be in the data. Gene sets with fewer genes overlapping with the dataset are dropped.
#' @param min_median_expression numeric, the minimum value for median expression of genes in a gene set. Any gene sets with median expression lower than this value are removed.
#' @param remove_missing_genes logical, whether to remove genes from the gene sets if they are not found in counts. Defaults to FALSE.
#' @export
#' @usage \code{filter_gene_sets(
#' gene_sets, counts,
#' counts_genes_by="rows",
#' min_genes=2,
#' min_median_expression=1,
#' remove_missing_genes=FALSE}
filter_gene_sets <-
  function(gene_sets, counts,
           counts_genes_by="rows",
           min_genes=2,
           min_median_expression=1,
           remove_missing_genes=FALSE) {
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
    
    # calculate # of genes present
    gene_sets.ngenes_present <-
      sapply(gene_sets, function(x) sum(x %in% rownames(counts)))
    gene_sets.median_exp <-
      sapply(gene_sets, function(x) median(counts[na.omit(match(x, rownames(counts))),]))
    
    # filter gene sets
    gene_sets <-
      gene_sets[gene_sets.ngenes_present >= min_genes &
                  gene_sets.median_exp >= min_median_expression]
    
    # convert back to original class, and fill in matrix or data.frame with NAs
    if (class.orig %in% c("matrix", "data.frame")) {
      seq_max <- seq_len(max(sapply(gene_sets, length)))
      gene_sets <- sapply(gene_sets, "[", i=seq_max)
      if (class.orig=="data.frame")
        gene_sets <- data.frame(gene_sets)
    }
  
    return(gene_sets)
  }
