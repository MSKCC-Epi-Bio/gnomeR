#' oncoKB.annotate
#'
#' OncoKB annotaor for MAF files.
#'
#' @param maf
#' @return OncoKb annotate MAF.
#' @export
#'
#' @examples

oncoKB.annotate <- function(input,output,lib.path = NULL){

  if(is.null(lib.path)) lib.path <- .libPaths()

  system(paste0("python ",lib.path,"/gnomeR/MafAnnotator.py -i ",input," -o ",output ))

}

