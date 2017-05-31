#' Calculate Fisher's exact test statistics for a gene set in a subset of a gene list
#'
#' Calculate Fisher's exact test statistics for a gene set in a subset of a list of genes, versus the
#' background to that subset (the full list of genes). i.e. test if the gene set is enriched in some
#' subset of a full gene list. The subset could be a WGCNA modules, a set of genes with significant
#' differential expression, etc. This tests against the null hypothesis of even distribution of the
#' contingency table; the test statistic is compared to a hypergeometric distribution. By default, it
#' tests for positive enrichment (genes in the gene set found more often in the subset than expected).
#' @param gene_set a character vector, the names of the genes in the gene set.
#' @param gene_list a character vector, the full gene list, including the subset to test.
#' @param index a logical or numeric vector. Either a logical vector with same length as \code{gene_list}, with TRUE/FALSE for each gene (where TRUE indicates membership in the subset), or a numeric vector with the indices from \code{gene_list} of the genes in the subset.
#' @param alternative the direction to test for enrichment. Passed to \code{fisher.test}. The default value, "greater", tests for members of \code{gene_set} found more often than expected in the specified subset of \code{gene_list}.
#' @param ... optional arguments to be passed to \code{fisher.test}
#' @export
#' @return A data frame with a single row, containing the p-value, the upper and lower confidence intervals, and the estimated odds ratio.
#' @usage \code{geneSetExactTest(gene_set, gene_list, index, alternative="greater", ...)}
geneSetExactTest <- function(gene_set, gene_list, index, alternative="greater", ...) {
  if (is.logical(index)) index <- which(index)
  
  cont.table <- matrix(nrow=2, ncol=2)
  cont.table[1,1] <- length(setdiff(gene_list[-index], gene_set))
  cont.table[1,2] <- length(intersect(gene_list[-index], gene_set))
  cont.table[2,1] <- length(setdiff(gene_list[index], gene_set))
  cont.table[2,2] <- length(intersect(gene_list[index], gene_set))
  
  test <- fisher.test(x=cont.table, alternative=alternative, ...)
  return(
    data.frame(
      p.value = test$p.value,
      conf.lower = test$conf.int[1],
      conf.upper = test$conf.int[2],
      odds_rat.estimate = test$estimate))
}
