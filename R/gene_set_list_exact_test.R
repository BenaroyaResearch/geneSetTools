#' Calculate Fisher's exact test statistics for a number of gene sets in a subset of a gene list
#'
#' Calculate Fisher's exact test statistics for a list of gene sets in a subset of a list of genes,
#' versus the background to that subset (the full list of genes). i.e. test if each gene set is enriched
#' in some subset of a full gene list. It is a wrapper for gene_set_exact_test, with multiple gene sets.
#' The subset could be a WGCNA modules, a set of genes with significant
#' differential expression, etc. This tests against the null hypothesis of even distribution of the
#' contingency table; the test statistic is compared to a hypergeometric distribution. By default, it
#' tests for positive enrichment (genes in the gene set found more often in the subset than expected).
#' @include gene_set_exact_test.R
#' @param gene_sets a list of character vectors, each vector containing the names of the genes in a single gene set.
#' @param gene_list a character vector, the full gene list, including the subset to test.
#' @param index a logical or numeric vector. Either a logical vector with same length as \code{gene_list}, with TRUE/FALSE for each gene (where TRUE indicates membership in the subset), or a numeric vector with the indices from \code{gene_list} of the genes in the subset.
#' @param alternative the direction to test for enrichment. Passed to \code{fisher.test} via \code{gene_set_exact_test}. The default value, "greater", tests for members of each gene set found more often than expected in the specified subset of \code{gene_list}.
#' @param fdr.method the FDR calculation method, passed to \code{p.adjust}. The default, "BH", specifies the Benjamini-Hochberg approach.
#' @param ... optional arguments to be passed to \code{fisher.test} via \code{gene_set_exact_test}.
#' @export
#' @return A data frame with a row for each gene set, containing the name of the gene set, the p-value, the upper and lower confidence intervals, and the estimated odds ratio.
#' @usage \code{gene_set_list_exact_test(gene_sets, gene_list, index,
#'   alternative="greater", fdr.method="BH", ...)}
gene_set_list_exact_test <-
  function(gene_sets, gene_list, index,
           alternative="greater", fdr.method="BH",
           sort_var="p.value", ...) {
    if (!(sort_var %in%
          c("gene_set", "p.value", "adj.p.value", "conf.lower", "conf.upper", "odds_rat.estimate")) &
        !(sort_var %in% 1:6))
      stop("sort_var must be a valid name or number of a column in results.")
    results <- data.frame(gene_set=character(), p.value=numeric(),
                          conf.lower=numeric(), conf.upper=numeric(),
                          odds_rat.estimate=numeric())
    for (i in 1:length(gene_sets)) {
      results[i,1] <- names(gene_sets)[i]
      results[i,2:5] <- gene_set_exact_test(gene_sets[[i]], gene_list, index, alternative=alternative, ...)
    }
    results$adj.p.value <- p.adjust(results$p.value, method=fdr.method)
    results <- results[,c("gene_set", "p.value", "adj.p.value", "conf.lower", "conf.upper", "odds_rat.estimate")]
    results <- results[order(results[,sort_var], decreasing=FALSE),]
    return(results)
  }
