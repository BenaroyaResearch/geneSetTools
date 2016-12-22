#' Reformat output of PCGSE for easier viewing
#' 
#' This function reformats the output of PCGSE, allowing for easier viewing of the results. It calculates the
#' adjusted p-values, and converts the list into a data frame with the gene set, test statistics, p-values, and
#' adjusted p-values as columns.
#' @param result the output from PCGSE. Usually a list with two matrics: p.values and statistics.
#' @param sort_by_PC the number of the PC to sort results by
#' @param gene_sets the gene sets with which PCGSE was run. Used for pulling in the gene set names if they are missing
#' @export
#' @usage \code{pcgse_process(result, sort_by_PC=1, gene_sets=NULL)}
#' @return A data frame, with one row for each gene set, and 3 columns for each PC tested by PCGSE. Columns include the test statistic, the p-value, and the adjsuted p-value.
pcgse_process <- function(result, sort_by_PC=1, gene_sets=NULL) {
  if (!is.null(gene_sets)) {
    if (is.list(gene_sets)) {
      gene_set_names <- names(gene_sets)
    } else if (is.matrix(gene_sets)) {
      gene_set_names <- rownames(gene_sets)
    } else {
      stop("Gene set object type not recognized.")
    }
  } else if (!is.null(rownames(result[[1]]))) {
    gene_set_names <- rownames(result[[1]])
  } else {
    cat("Gene set names not found; setting all to 'NA'\n")
    gene_set_names <- rep(NA, nrow(result[[1]]))
  }
  
  result.processed <- data.frame(
    gene_set = gene_set_names
  )
  if (is.null(result$adj.p.values))
    result$adj.p.values <- matrix(p.adjust(result$p.values), ncol=ncol(result$p.values))
  for (i in 1:ncol(result[[1]])) {
    result.processed[[paste0("PC",i,".statistic")]] <- result$statistics[,i]
    result.processed[[paste0("PC",i,".p.values")]] <- result$p.values[,i]
    result.processed[[paste0("PC",i,".adj.p.values")]] <- result$adj.p.values[,i]
  }
  result.processed <- result.processed[order(result.processed[,paste0("PC",sort_by_PC,".p.values")]),]
}
