#' Read a gene set in .gmt file format
#' 
#' This function reads a gene set in .gmt file format into a list. By default, it uses the first
#' element of each gene set as the name, and assumes the first gene in each set is the third element.
#' @param filename path to the .gmt file
#' @param out_format the format to output, either "list" or "matrix".
#' @param empty_val the value to use for empty cells. Ignored if out_format == "list".
#' @param name_loc index (position) of the gene set names in each gene set. List elements are assigned these as names. Defaults to 1.
#' @param genes_start index (position) of the first gene name in each gene set. Allows skipping descriptions and extra elements. Defaults to 3.
#' @export
#' @usage \code{read_gmt(filename, out_format="matrix", empty_val=NA, name_loc=1, genes_start=3)}
read_gmt <- function(filename, out_format="matrix", empty_val=NA, name_loc=1, genes_start=3) {
  out_format==match.arg(out_format, choices=c("matrix", "list"))
  if (out_format=="list") {
    data <- readLines(filename)
    data <- strsplit(data, "\t")
    data <- lapply(data, FUN=na.omit) # remove NAs
    names(data) <- sapply(data, FUN=function(x) x[name_loc])
    data <- lapply(data, FUN=function(x) x[genes_start:length(x)])
  } else if (out_format=="matrix") {
    data <- read.delim(filename, header=FALSE, stringsAsFactors=FALSE)
    data <- data[stringr::str_detect(data[,2], "http"),]
    data <- t(data)
    colnames(data) <- data[1,]
    data <- data[-(1:2),]
    data[data==""] <- empty_val
  }
  return(data)
}
