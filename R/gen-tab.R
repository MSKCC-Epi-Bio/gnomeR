#' gen.tab
#'
#' Creates a summary table of the distribution of the genetic features by a specific outcome/covariate of interest. The difference
#' is tested using Fisher's exact test and further adjusted for multiple comparisons. Note that continuous genetic factors
#' are dichotomized at their median.
#'
#' @param gen.dat A matrix or dataframe, with patients as rows and features as columns.
#' @param outcome A leveled vector of length equal to the number of rows in gen.dat.
#' @return p : an oncoprint object
#'
#' @export
#'
#' @examples library(gnomeR)
#' patients <- unique(clin$DMPID)[1:500]
#' mut.only <- create.bin.matrix(patients = patients,maf = mut)
#' gen.dat <- mut.only$mut
#' outcome = as.factor(as.character(clin$Sex[match(patients,clin$DMPID)]))
#' test <- gen.tab(gen.dat,outcome)
#' head(test)

gen.tab <- function(gen.dat,outcome){
  fits <- as.data.frame(t(apply(gen.dat,2,function(x){
    #print(x)
    if(length(unique(x)) > 20){x <- ifelse(x>median(x,na.rm = T) ,1 ,0)}
    test <- fisher.test(x,outcome)
    out <- c()
    for(i in 1:length(levels(outcome))){
      out <- c(out,sum(x[which(outcome == levels(outcome)[i])]))
    }
    out <- c(out,test$estimate,test$p.value)
    if(!is.null(test$estimate)){
      names(out) <- c(levels(outcome)[1:length(levels(outcome))],"OddsRatio","Pvalue")}
    else{names(out) <- c(levels(outcome)[1:length(levels(outcome))],"Pvalue")}
    return(out)
  })))

  fits$FDR <- p.adjust(fits$Pvalue,method="fdr")
  fits <- fits[order(fits$Pvalue),]
  return(fits)
}
