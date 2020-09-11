#' custom_pathway
#' Enables creation of a custom pathway binary matrix from a binmat() `object`
#' @param mat a binmat() binary matrix
#' @param pathway a dataframe/tibble with first column the name of the genes, and second column their corresponding
#' pathway
#' @return pathway_dat : a binary matrix of pathway level alterations
#' @export
#' @examples library(gnomeR)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,maf = mut)
#' pathway <- as.data.frame(cbind(c("PIK3CA","KRAS","TERT","TP53"),c("path1","path1","path2","path3")))
#' custom_pathway(mat = bin.mut, pathway = pathway)
#' @import dplyr
#' @import dtplyr
#' @import stringr


custom_pathway <- function(mat, pathway){

  colnames(pathway) <- c("hugo_symbol","pathway")
  pathway <- pathway %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol),
           pathway = as.character(.data$pathway))

  pathway_dat <- as.data.frame(do.call('cbind',lapply(unique(pathway$pathway[!is.na(pathway$pathway)]),function(x){
    genes <- as.character(unlist(pathway %>%
                                   filter(.data$pathway == x) %>% select(.data$hugo_symbol)))
    as.numeric(apply(mat %>% select(starts_with(genes)),1,function(y){
      ifelse(sum(abs(as.numeric(as.character(y))),na.rm = T)>0, 1,0)
    }))
  })))
  colnames(pathway_dat) <- unique(pathway$pathway[!is.na(pathway$pathway)])
  rownames(pathway_dat) <- rownames(mat)
  return(pathway_dat)
}

