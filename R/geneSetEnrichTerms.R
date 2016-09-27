#' Extract enrichment terms from a list of gene sets significantly up and/or down in a gene set analysis
#'
#' This is a fetch function; it generates vectors of terms associated with gene sets that are
#' significantly up or down in a gene set analysis. Basically, it takes the output of the gene set
#' analysis, identifies gene sets that are significantly up and/or down in the that analysis, and
#' extracts terms related to those gene sets from a data frame describing the gene sets.
#' @param results the output from a gene set analysis (e.g. roast, mroast, camera, gene)
#' @param gene_set_df a data frame containing the details of the gene set annotations.
#' @param method a string, the method used to run the gene set analysis. This provides the function with the expected structure of \code{results}, so that various other parameters do not need to be specified. Current options include "camera", "mroast", "geneSetListExactTest", and "WGCNAmodulesgeneSetListExactTest
#' @param results.name_col name or number of the column of \code{results} containing the gene set names. If gene set names are the row names, specify "rownames".
#' @param results.sig_col name or number of the column of \code{results} containing the gene set significance numbers.
#' @param results.direction_col name number of the column of \code{results} containing the direction (Up/Down).
#' @param gene_set_df.name_col name or number of the column of \code{gene_set_df} containing the gene set names. If gene set names are the row names, specify "rownames".
#' @param gene_set_df.terms_col name or number of the column of \code{gene_set_df} containing the gene set terms or descriptions.
#' @param threshold numeric, the maximum value of the significance statistic for which to include gene sets. Defaults to 0.01.
#' @param direction character, the direction of gene sets to include. Can be "up", "down", or "both". Only used if results object contains a column for direction.
#' @import stringr
#' @export
#' @return A data frame (or list of data frames, one for each direction) containing the names of the enriched gene sets and the terms or description associated with each gene set.
#' @usage \code{
#' geneSetEnrichTerms(results, gene_set_df, method=NULL, 
#'                    results.name_col=NULL, results.sig_col=NULL, results.direction_col=NULL,
#'                    gene_set_df.name_col, gene_set_df.terms_col,
#'                    threshold=0.01, direction=NULL)}
geneSetEnrichTerms <- function(results, gene_set_df, method=NULL, 
                               results.name_col=NULL, results.sig_col=NULL, results.direction_col=NULL,
                               gene_set_df.name_col, gene_set_df.terms_col,
                               threshold=0.01, direction=NULL) {
  if (is.null(method) &
      (is.null(results.name_col) | is.null(results.sig_col) | is.null(results.direction_col)))
    stop("Must provide either 'method' or all column specifiers for 'results'")
  
  if (method %in% c("camera", "mroast")) {
    if (is.null(results.name_col)) results.name_col <- "rownames"
    if (is.null(results.sig_col)) results.sig_col <- "FDR"
    if (is.null(results.direction_col)) results.direction_col <- "Direction"
  } else if (method == "geneSetListExactTest") {
    if (is.null(results.name_col)) results.name_col <- "gene_set"
    if (is.null(results.sig_col)) results.sig_col <- "adj.p.value"
  } else if (method == "WGCNAmodulesGeneSetListExactTest") {
    # rerun geneSetEnrichTerms for each module, and return a list with results for each module
    output <- list()
    for (i in unique(results[["module"]])) {
      output[[i]] <-
        geneSetEnrichTerms(results[results[["module"]]==i,],
                           gene_set_df=gene_set_df, method="geneSetListExactTest",
                           gene_set_df.name_col=gene_set_df.name_col,
                           gene_set_df.terms_col=gene_set_df.terms_col,
                           results.name_col=results.name_col,
                           results.sig_col=results.sig_col,
                           threshold=threshold)
    }
    return(output)
  }
  
  if (results.name_col %in% c("rownames", "row.names"))
    results[,results.name_col] <- rownames(results)
  if (gene_set_df.name_col %in% c("rownames", "row.names"))
    gene_set_df[,gene_set_df.name_col] <- rownames(gene_set_df)
  
  if (is.null(direction)) {
    gene_sets.tmp <- results[results[[results.sig_col]] <= threshold,
                             results.name_col]
    terms <- data.frame(
      gene_set = gene_sets.tmp,
      terms = gene_set_df[match(gene_sets.tmp, gene_set_df[[gene_set_df.name_col]]),
                          gene_set_df.terms_col])
  } else {
    direction <- str_to_lower(direction)
    results[,results.direction_col] <- sapply(results[,results.direction_col], str_to_lower)
    
    if (direction=="both") direction <- c("up","down")
    
    terms <- list()
    for (i in direction) {
      gene_sets.tmp <- results[(results[[results.direction_col]]==i) &
                                 (results[[results.sig_col]] <= threshold),
                               results.name_col]
      terms[[i]] <- data.frame(
        gene_set = gene_sets.tmp,
        terms = gene_set_df[match(gene_sets.tmp, gene_set_df[[gene_set_df.name_col]]),
                            gene_set_df.terms_col])
    }
  }
  
  return(terms)
}