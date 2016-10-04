#' Generate a binary incidence matrix for gene sets and counts
#' 
#' This function creates a binary or index-based gene set incidence matrix. It expects gene sets as named lists of character vectors.
#' @param gene_sets the gene sets. Either a list of gene sets, with each as a character vector, or a matrix, with one gene set per column, and the gene set names as column names. Names of list elements are used as gene set names.
#' @param counts a matrix or data frame of counts. The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}.
#' @param counts_genes_by the dimension of \code{counts} corresponding to genes. Can be "rows" or "columns" or partial matches.
#' @param gene_sets_by the dimension of the result that should correspond to gene sets. Can be "rows" or "columns" or partial matches. Not used if method is "index".
#' @param min_genes an integer, the minimum number of genes from a gene set that must be in the data. Gene sets with fewer genes overlapping with the dataset are dropped.
#' @param method character, the type of incidence object to return. Either "binary" for a binary matrix or "index" for a list of vector-based indices.
#' @export
#' @usage \code{gene_sets_incMatrix(
#' gene_sets, counts,
#' counts_genes_by="rows",
#' gene_sets_by="rows",
#' min_genes=2,
#' method="binary"}
gene_sets_incMatrix <- 
  function(gene_sets, counts,
           counts_genes_by="rows",
           gene_sets_by="rows",
           min_genes=2,
           method="binary") {
  counts_genes_by <- match.arg(counts_genes_by, choices=c("rows", "columns"))
  gene_sets_by <- match.arg(gene_sets_by, choices=c("rows", "columns"))
  method <- match.arg(method, choices=c("index", "binary"))
  if (counts_genes_by=="columns") counts <- t(counts)
  
  if (is.data.frame(gene_sets) | is.matrix(gene_sets)) {
    if (gene_sets_by=="columns") gene_sets <- t(gene_sets)
    if (method=="index") {
      incMatrix <- list()
      for (i in colnames(gene_sets))
        incMatrix[[i]] <- na.omit(match(gene_sets[,i]), rownames(counts))
    } else if (method=="binary") {
      incMatrix <- t(apply(gene_sets, 1, function(x) {rownames(counts) %in% x}))
      colnames(incMatrix) <- rownames(counts)
      incMatrix <- incMatrix + 0
      incMatrix <- incMatrix[rowSums(incMatrix) >= min_genes,]
    }
    
  } else if (is.list(gene_sets)) {
  
    if (method=="index") {
      incMatrix <- list()
      for (i in names(gene_sets))
        incMatrix[[i]] <- na.omit(match(gene_sets[[i]], rownames(counts)))
    } else if (method=="binary") {
      incMatrix <- matrix(data=as.logical(NA), nrow=length(gene_sets), ncol=nrow(counts),
                          dimnames=list(names(gene_sets), rownames(counts)))
      n_gene_sets <- length(gene_sets)
      for (i in 1:n_gene_sets) {
        incMatrix[i,] <- rownames(counts) %in% gene_sets[[i]]
      }
      incMatrix <- incMatrix + 0
      incMatrix <- incMatrix[rowSums(incMatrix) >= min_genes,]
      
      if (gene_sets_by=="columns") incMatrix <- t(incMatrix)
    }
  }
    
  return(incMatrix)
}
