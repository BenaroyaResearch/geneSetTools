#' Write out counts to a .txt file compatible with GSEA
#'
#' Write out expression counts to a .txt file, for use with the Java implementation of GSEA. This function
#' formats the output to match the .txt counts file format used by GSEA.
#' @param filename A string, the name of the file to output
#' @param counts A matrix or data frame, the counts to be used in GSEA, with samples in columns and genes in rows. Should have rows labeled with gene identifiers.
#' @export
#' @usage \code{write_GSEA_counts(filename, counts)}
write_GSEA_counts <- function(filename, counts) {
  data <- data.frame(NAME=rownames(counts),
                     DESCRIPTION=rownames(counts),
                     counts)
  write.table(data, file=filename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}