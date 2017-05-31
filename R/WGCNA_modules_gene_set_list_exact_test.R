#' Calculate Fisher's exact test statistics for a number of gene sets in WGCNA modules
#'
#' Calculate Fisher's exact test statistics for a list of gene sets in WGCNA modules, versus the 
#' background of the full list of genes that the WGCNA modules were constructed from. i.e. test if 
#' each gene set is enriched in each module. It is a wrapper for gene_set_list_exact_test, with multiple 
#' subsets of the full gene list. This tests against the null hypothesis of even distribution of the
#' contingency table; the test statistic is compared to a hypergeometric distribution. By default, it
#' tests for positive enrichment (genes in the gene set found more often in the subset than expected).
#' @include gene_set_list_exact_test.R
#' @param gene_sets a list of character vectors, each vector containing the names of the genes in a single gene set.
#' @param gene_modules a data frame containing gene names ("gene"), WGCNA modules ("color_assigned"), and connectivity to the assigned module ("kME_color_assigned"). It should contain all genes included in the module construction, as the full set in \code{gene_modules} is used as the background.
#' @param threshold a scalar, the value of kME to use as a threshold for inclusion of genes in WGCNA modules. Genes with connectivity below this value are considered to NOT be members of the modules; they are included in the background list. Based on some minimal testing, threshold does not seem to have a huge effect on which gene sets are near the top of the list, but it does affect how many gene sets are found significant.
#' @param alternative the direction to test for enrichment. Passed to \code{fisher.test}. The default value, "greater", tests for members of each gene set found more often than expected in each WGCNA module.
#' @param fdr.method the FDR calculation method, passed to \code{p.adjust}. The default, "BH", specifies the Benjamini-Hochberg approach.
#' @param p.adjust_set the set of p-values to include in FDR q-value (adjusted p-value) calculation. The default, "all", aggregates all tests of all gene sets against all WGCNA modules, and calculates adjusted p-values for that set. Other possible values: "module" adjusts p-values for all gene sets tested within each module. Other values result in a warning, and adjustment as for "module".
#' @param return_format the format in which to return the results. The default, "df", returns a data frame, with a column for module. "list" returns a list, with each element containing the results for a single module. Other values result in an error.
#' @param return_sort the name(s) of the column(s) to sort results on. Defaults to "p.value".
#' @param ... optional arguments to be passed to \code{fisher.test}.
#' @export
#' @return Either a data frame with a row for each combination of module and gene set; or a list of data frames, one for module, with each containing rows for the results of each gene set tested. In either case, the data frame(s) contain(s) the name of the gene set, the p-value, the upper and lower confidence intervals, and the estimated odds ratio.
#' @usage \code{WGCNA_modules_gene_set_list_exact_test(gene_sets, gene_modules, threshold=0.75,
#'   alternative="greater", fdr.method="BH", p.adjust_set="all",
#'   return_format="df", return_sort="p.value", ...)}
WGCNA_modules_gene_set_list_exact_test <-
  function(gene_sets, gene_modules, threshold=0.75,
           alternative="greater", fdr.method="BH", p.adjust_set="all",
           return_format="df", return_sort="p.value", ...) {
    results <- data.frame(gene_set=character(), p.value=numeric(),
                          adj.p.value=numeric(),
                          conf.lower=numeric(), conf.upper=numeric(),
                          odds_rat.estimate=numeric())
    
    for (i in unique(gene_modules$color_assigned)) {
      index.tmp <- (gene_modules$color_assigned==i &
                      gene_modules$kME_color_assigned >= threshold)
      results.tmp <-
        gene_set_list_exact_test(gene_sets, gene_list=gene_modules$gene,
                             index=index.tmp, ...)
      results <- rbind(results, cbind(rep(i,nrow(results.tmp)), results.tmp))
    }
    names(results)[1] <- "module"
    
    if (p.adjust_set=="all") {
      results$adj.p.value <- p.adjust(results$p.value, method=fdr.method)
    } else if (p.adjust_set!="module") {
      cat("Warning: p-value adjustment set not recognized. Adjusting by module, but not across modules.")
    }
    
    results <- results[order(results[,return_sort]),]
    
    # return results in df or list format
    if (return_format=="df") {
      return(results) 
    } else if (return_format=="list") {
      results.list <- list()
      for (i in unique(results$module)) {
        results.list[[i]] <- results[results$module==i,]
      }
      return(results.list)
    } else stop("Error: return format not recognized. Please specify 'df' or 'list'.")
  }
