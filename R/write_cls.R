#' Write out phenotypes to a .cls file compatible with GSEA
#'
#' Write out phenotype values to a .cls file, for use with the Java implementation of GSEA. Samples are
#' output in the order they are found in an associated counts object so that they are properly associated
#' in GSEA.
#' @param file_prefix A string, the (optional) path and file name to output. The name of the variable, and ".cls" will be appended to this.
#' @param counts A matrix or data frame, the counts to be used in GSEA, with samples in columns and genes in rows. Should have columns labeled with sample identifiers. Necessary in order to determine the proper order of sample phenotype values in the output file.
#' @param design A data frame containing (at minimum) columns for the sample identifiers and phenotype values.
#' @param var_to_test A string, the name of the column in \code{design} that contains the phenotype values. Also added to the output file name.
#' @param libID_col Optional, the name or number of the column in \code{design} that contains sample identifiers matching the column names in \code{counts}. If "row.names" is specified, the rownames of \code{design} are used.
#' @export
#' @usage \code{write_cls(file_prefix, counts, design, var_to_test, libID_col="lib.id")}
write_cls <- function(file_prefix, counts, design, var_to_test, libID_col="lib.id") {
  if (is.matrix(design)) design <- as.data.frame(design)
  if (libID_col == "row.names" & !("row.names" %in% colnames(design))) design$row.names <- rownames(design)
  if (!(libID_col %in% colnames(design))) stop(paste0("Design object is missing column ", libID_col, ", where I expected to find library identifiers."))
  if (!setequal(colnames(counts), design[,libID_col])) stop("Library identifiers in the counts and design objects do not match. Check that the objects contain the same libraries.")
  if (any(duplicated(colnames(counts))) | any(duplicated(design[,libID_col]))) stop("Duplicated library identifiers in counts or design object.")
  
  sink(paste(file_prefix, var_to_test, "cls", sep="."))
  if (is.numeric(design[,var_to_test])) {
    cat("#numeric\n")
    cat("#", var_to_test, "\n", sep="")
    cat(design[match(colnames(counts), design[,libID_col]), var_to_test]); cat("\n")
  } else {
    design[,var_to_test] <- factor(design[,var_to_test])
    cat(ncol(counts), length(levels(design[,var_to_test])), "1\n", sep=" ")
    cat("#", levels(design[,var_to_test]), sep=" "); cat("\n")
    cat(as.numeric(design[match(colnames(counts), design[,libID_col]), var_to_test])-1); cat("\n")
  }
  sink()
}