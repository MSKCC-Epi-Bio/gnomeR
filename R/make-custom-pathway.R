#' custom_pathway
#' \%lifecycle{stable}
#' Enables creation of a custom pathway binary matrix from a binmat() `object`. Similarly to the internal file 'pathways.csv', this function takes as input
#' a data frame containing the name of the pathways of interest, and for each pathways a character vector of the genes of interest. Note that the different
#' events to be considered in each pathways must be considered separetely. For example, if one wishes to consider TP53 deletions in a given pathway, one must
#' specify "TP53.Del" in the character vector for that pathway.
#' @param mat a binmat() binary matrix
#' @param pathway a dataframe/tibble with first column the name of the genes, and second column their corresponding
#' pathway
#' @return pathway_dat : a binary matrix of pathway level alterations
#' @export
#' @examples library(gnomeR)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,maf = mut)
#' pathway <- as.data.frame(cbind(c("path1","path1","path2","path3"),c("PIK3CA","KRAS","TERT","TP53")))
#' custom_pathway(mat = bin.mut, pathway = pathway)
#' ### Considering CNA as well ###
#' bin.mut <- binmat(patients = patients,maf = mut,cna = cna)
#' pathway <- as.data.frame(cbind(c("path1","path1","path2","path3","path3"),c("PIK3CA","KRAS","TERT","TP53","TP53.Del")))
#' custom_pathway(mat = bin.mut, pathway = pathway)
#' @import dplyr
#' @import dtplyr
#' @import stringr


custom_pathway <- function(mat, pathway){

  colnames(pathway) <- c("pathway","hugo_symbol")
  pathway <- pathway %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol),
           pathway = as.character(.data$pathway))

  pathway_dat <- as.data.frame(
    do.call('cbind',
            lapply(unique(pathway$pathway[!is.na(pathway$pathway)]),function(x){
              genes <- unlist(strsplit(as.character(unlist(pathway %>%
                                                      filter(.data$pathway == x) %>% select(.data$hugo_symbol))),","))#[[1]]
              if(all(is.na(match(genes,colnames(mat))))){
                warning("None of the names specified in pathway: ",x," were found.")
                return(rep(0, rnow(mat)))
              }

              as.numeric(apply(
                as.data.frame(mat[,na.omit(match(genes,colnames(mat)))]),1,function(y){ #%>% select(starts_with(genes))
                ifelse(sum(abs(as.numeric(as.character(y))),na.rm = T)>0, 1,0)
              }))
            })))
  colnames(pathway_dat) <- unique(pathway$pathway[!is.na(pathway$pathway)])
  rownames(pathway_dat) <- rownames(mat)
  return(pathway_dat)
}

