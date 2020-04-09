#' oncoKB.annotate
#'
#' OncoKB annotator for MAF files.
#'
#' @param input A MAF file to be annotated.
#' @param output Name of the output file
#' @param lib.path Path to the
#' @return OncoKb annotate MAF.
#' @export
#'
#' @examples
#'
#' # Needs to be run on LINUX with python
#' # oncoKB.annotate(input = "ExampleMaf", out = "ExampleOut")

oncoKB.annotate <- function(input,output,lib.path = NULL){

  if(is.null(lib.path)) lib.path <- .libPaths()

  system(paste0("python ",lib.path,"/gnomeR/MafAnnotator.py -i ",input," -o ",output ))

}

