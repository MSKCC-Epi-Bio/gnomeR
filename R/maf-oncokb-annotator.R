#' oncoKB.annotate
#'
#' OncoKB annotaor for MAF files.
#'
#' @param maf
#' @return OncoKb annotate MAF.
#' @export
#'
#' @examples
#'
#' @import

oncoKB.annotate <- function(input,output,lib.path = NULL){

  if(is.null(lib.path)) lib.path <- .libPaths()

  system(paste0("python ",lib.path,"/python/MafAnnotator.py -i ",input," -o ",output ))

}

