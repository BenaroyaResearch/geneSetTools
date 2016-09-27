#' Remove duplicate elements from a list of vectors
#'
#' Remove duplicate elements from a list of vectors. This cycles through a list of vectors, repeatedly
#' checking the elements of each vector against all other vectors in the list. If any elements are
#' shared with other vectors, one is removed (can be the first or last such duplicate in the vector),
#' and the algorithm moves to the next vector. This continues until no elements are included in more than
#' one of the vectors in the list. It is used in the context of gene set analyses to remove genes that
#' are duplicated across multiple gene sets, so that the gene sets are non-overlapping.
#' @param data a list of character vectors.
#' @param remove which duplicate to remove; "first" always removes the first element found in another vector, "last" removes the last, and "random" selects one at random to remove
#' @export
#' @return A data frame with a row for each gene set, containing the name of the gene set, the p-value, the upper and lower confidence intervals, and the estimated odds ratio.
#' @usage \code{deduplicate_list_elements(data, remove="last")}
deduplicate_list_elements <- function(data, remove="last") {
  if (!(remove %in% c("last", "first", "random"))) stop("Removal order not recognized. Please specify 'first', 'last', or 'random'.")
  while (any(duplicated(unlist(data)))) {
    for (i in 1:length(data)) {
      duplicates.tmp <- data[[i]] %in% unlist(data[-i])
      
      if (any(duplicates.tmp)) {
        data[[i]] <- switch(remove,
                            first = data[[i]][-which(duplicates.tmp)[1]],
                            last = data[[i]][-which(duplicates.tmp)[sum(duplicates.tmp)]],
                            random = data[[i]][-sample(which(duplicates.tmp), size=1)],
                            data[[i]])
      }
    }
  }
  return(data)
}
  