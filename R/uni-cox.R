#' uni.cox
#'
#' Performs univariate cox proportional hazard model on every feature
#'
#' @param X Matrix/surv.datframe of genomic features, continuous or binary (note cannot handle categorical surv.dat for the moment).
#' @param surv.dat a surv.dat frame containing the survival information. This can be made of 2 or 3 columns. 1 or 2 for time,
#' and one for status (where 1 is event and 0 is no event).
#' @param surv.formula a survival formula with names matching those in surv.dat eg: Surv(time,status)~.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of patients
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @return tab A table of all the fits performed sorted by adjusted pvalues.
#' @return p An interactive plot of log(pvalue) by hazard ration.
#' @export
#'
#' @examples library(gnomeR)
#' patients <- unique(clin$DMPID)[1:100]
#' mut.only <- create.bin.matrix(patients = patients,maf = mut)
#' gen.dat <- mut.only$mut
#' surv.dat <-clin[match(patients,clin$DMPID),match(c("time","status") ,colnames(clin))]
#' surv.dat$status <- ifelse(surv.dat$status == "DECEASED",1,0)
#' surv.dat$time <- as.numeric(as.character(surv.dat$time))
#'
#' cox.fits <- uni.cox(X = gen.dat,surv.dat = surv.dat, surv.formula = Surv(time,status)~.)
#'
#' @import
#' dplyr
#' plotly
#' survival


uni.cox <- function(X,surv.dat,surv.formula,filter = 0){

  # filtering #
  if(!(filter >= 0 && filter < 1))
    stop("Please select a filter value between 0 and 1")
  if(filter > 0){
    # get binary cases #
    temp <- apply(X, 2, function(x){length(unique(x)) == 2})
    genes.bin <- names(temp[which(temp)])
    rm <- apply(X[,genes.bin], 2, function(x){sum(x)/length(x) < filter})
    genes.rm <- names(rm[which(rm)])
    X <- X %>%
      select(-one_of(genes.rm))
  }


  # appropriate formula
  survFormula <- as.formula(surv.formula)
  survResponse <- survFormula[[2]]

  if(!length(as.list(survResponse)) %in% c(3,4)){
    stop("ERROR : Response must be a 'survival' object with 'Surv(time, status)' or 'Surv(time1, time2, status)'.")
  }
  ### reprocess surv.dat
  if(length(as.list(survResponse)) == 3){
    colnames(surv.dat)[match(as.list(survResponse)[2:3],colnames(surv.dat))] <- c("time","status")
    LT = FALSE
    timevars <- match(c("time","status"),colnames(surv.dat))
    surv.dat <- surv.dat[,c(timevars,
                            c(1:ncol(surv.dat))[-timevars])]
  }
  if(length(as.list(survResponse)) == 4){
    colnames(surv.dat)[match(as.list(survResponse)[2:4],colnames(surv.dat))] <- c("time1","time2","status")
    LT = TRUE
    timevars <- match(c("time1","time2","status"),colnames(surv.dat))
    surv.dat <- surv.dat[,c(timevars,
                            c(1:ncol(surv.dat))[-timevars])]
  }

  #################################



  ##### univariate volcano plot #####
  if(LT == F){
    uni <- as.data.frame(t(apply(X,2,function(x){
      fit <- coxph(Surv(surv.dat$time,surv.dat$status)~x)
      out <- as.numeric(summary(fit)$coefficients[c(1,5)])
      return(out)
    })))
  }

  if(LT == T){
    uni <- as.data.frame(t(apply(X,2,function(x){
      fit <- coxph(Surv(surv.dat$time1,surv.dat$time2,surv.dat$status)~x)
      out <- as.numeric(summary(fit)$coefficients[c(1,5)])
      return(out)
    })))
  }
  colnames(uni) <- c("Coefficient","Pvalue")
  uni$Feature <- colnames(X)
  uni$FDR <- p.adjust(uni$Pvalue, method = "fdr")
  uni <- uni[order(uni$Pvalue),]

  uniVolcano <- plot_ly(data = uni, x = ~Coefficient, y = ~-log10(Pvalue),
                        text = ~paste('Feature :',Feature,
                                      '<br> Hazard Ratio :',round(exp(Coefficient),digits=3)),
                        mode = "markers") %>%
    layout(title ="Volcano Plot")


  return(list("tab" = uni,"p"=uniVolcano))
}

