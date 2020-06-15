#' uni.cox
#' Performs univariate cox proportional hazard model on every feature
#' @param X Matrix/surv.datframe of genomic features, continuous or binary (note cannot handle categorical surv.dat for the moment).
#' @param surv.dat a surv.dat frame containing the survival information. This can be made of 2 or 3 columns. 1 or 2 for time,
#' and one for status (where 1 is event and 0 is no event).
#' @param surv.formula a survival formula with names matching those in surv.dat eg: Surv(time,status)~.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of patients
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @param genes a character vector of gene names that will be the only ones to be kept. Default is NULL, all genes are used.
#' @return tab A table of all the fits performed sorted by adjusted pvalues.
#' @return p An interactive plot of log(pvalue) by hazard ration.
#' @return KM List of survival plots of the top 10 most significant genes
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' gen.dat <- binmat(patients = patients,maf = mut)
#' surv.dat <- clin.patients %>%
#' filter(X.Patient.Identifier %in%
#' abbreviate(patients,strict = TRUE, minlength = 9)) %>%
#'   select(X.Patient.Identifier,Overall.Survival..Months.,
#'    Overall.Survival.Status) %>%
#'   rename(DMPID = X.Patient.Identifier,
#'    time = Overall.Survival..Months.,
#'    status = Overall.Survival.Status) %>%
#'   mutate(time = as.numeric(as.character(time)),
#'          status = ifelse(status == "LIVING",0,1)) %>%
#'   filter(!is.na(time))
#' X <- gen.dat[match(surv.dat$DMPID,
#' abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9)),]
#' uni.cox(X = X, surv.dat = surv.dat,
#' surv.formula = Surv(time,status)~.,filter = 0.05)
#' @import
#' dplyr
#' survival
#' survminer
#' @importFrom plotly plot_ly layout


uni.cox <- function(X,surv.dat,surv.formula,filter = 0,genes = NULL){

  # filtering #
  if(!(filter >= 0 && filter < 1))
    stop("Please select a filter value between 0 and 1")
  if(!is.null(genes) && sum(colnames(X) %in% genes) == 0)
    stop("The genes argument inputted did not match any of the columns in the features matrix X.")
  else if(!is.null(genes) && sum(colnames(X) %in% genes) > 0){
    genes <- genes[genes %in% colnames(X)]
    X <- as.data.frame(X %>%
                         select(genes))
  }

  if(is.null(dim(X)) )
    stop("Only one or fewer genes were found from the 'genes' argument. We need a minimum of two.")
  # remove all columns that are constant #
  if(length(which(apply(X, 2, function(x){length(unique(x[!is.na(x)])) == 1} || all(is.na(x))))) > 0)
    X <- X[,-which(apply(X, 2, function(x){length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))}))]

  if(filter > 0){
    rm <- apply(X, 2, function(x) {
      # if(is.numeric(x))
      #   sum(x, na.rm = T)/length(x) < filter
      # else
      any(summary(as.factor(x[!is.na(x)]))/length(x[!is.na(x)]) < filter)
    })
    genes.rm <- names(rm[which(rm)])
    length(genes.rm)
    X <- X %>% select(-one_of(genes.rm))
  }
  if(is.null(dim(X)) )
    stop("Only one or fewer genes are left after filtering. We need a minimum of two. Please relax the filter argument.")

  # appropriate formula
  survFormula <- stats::as.formula(surv.formula)
  survResponse <- survFormula[[2]]

  if(!length(as.list(survResponse)) %in% c(3,4)){
    stop("ERROR : Response must be a 'survival' object with 'Surv(time, status)~.' or 'Surv(time1, time2, status)~.'.")
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

  ###################################
  ##### univariate volcano plot #####
  if(!LT){
    uni <- lapply(X,function(x){
      if(length(unique(x[!is.na(x)]))==2){
        x <- as.numeric(as.character(x))
        fit <- coxph(Surv(surv.dat$time,surv.dat$status)~x)
        out <- c(as.numeric(summary(fit)$coefficients[c(1,5)]),sum(x)/length(x))
        names(out) <- c("Coefficient","Pvalue","MutationFrequency")
      }
      else{
        x <- factor(as.numeric(as.character(x)),
                    levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(x)))])
        fit <- coxph(Surv(surv.dat$time,surv.dat$status)~x)
        out <- cbind(summary(fit)$coefficients[,c(1,5)],(summary(x)/length(x))[-1])
        rownames(out) <- levels(x)[-1]
        colnames(out) <- c("Coefficient","Pvalue","MutationFrequency")
      }
      return(out)
    })
  }

  if(LT){
    uni <- lapply(X,function(x){
      if(length(unique(x[!is.na(x)]))==2){
        x <- as.numeric(as.character(x))
        fit <- coxph(Surv(surv.dat$time1,surv.dat$time2,surv.dat$status)~x)
        out <- c(as.numeric(summary(fit)$coefficients[c(1,5)]),sum(x)/length(x))
        names(out) <- c("Coefficient","Pvalue","MutationFrequency")
      }
      else{
        x <- factor(as.numeric(as.character(x)),
                    levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(x)))])
        fit <- coxph(Surv(surv.dat$time1,surv.dat$time2,surv.dat$status)~x)
        out <- cbind(summary(fit)$coefficients[,c(1,5)],(summary(x)/length(x))[-1])
        rownames(out) <- levels(x)[-1]
        colnames(out) <- c("Coefficient","Pvalue","MutationFrequency")
      }
      return(out)
    })
  }

  ps <- c()
  for(i in 1:length(uni)){
    if(is.null(dim(uni[[i]])))
      ps[i] <- uni[[i]][2]
    else{
      ps[i] <- min(uni[[i]][,2])
      rownames(uni[[i]]) <- paste(names(uni)[[i]],rownames(uni[[i]]))
    }
  }
  uni <- do.call('rbind',uni[order(as.numeric(ps))])
  uni <- as.data.frame(cbind(rownames(uni),uni))
  colnames(uni) <- c("Feature","Coefficient","Pvalue","MutationFrequency")
  uni$FDR <- stats::p.adjust(as.numeric(as.character(uni$Pvalue)), method = "fdr")
  uni <- uni %>%
    mutate(
      Pvalue = as.numeric(formatC(stats::p.adjust(as.numeric(as.character(.data$Pvalue)),
                                       method = "fdr"), format = "e", digits = 2)),
      HR = as.numeric(exp(as.numeric(as.character(.data$Coefficient)))),
      Coefficient = as.numeric(formatC(as.numeric(as.character(.data$Coefficient)), format = "e", digits = 2)),
      HR = as.numeric(formatC(as.numeric(as.character(.data$HR)), format = "e", digits = 2)),
      MutationFrequency = as.numeric(round(as.numeric(as.character(.data$MutationFrequency)), digits = 4))
    )
  uniVolcano <- plot_ly(data = uni, x = ~Coefficient, y = ~-log10(Pvalue),
                        text = ~paste('Feature :',Feature,
                                      '<br> Hazard Ratio :',HR),
                        mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>%
    layout(title ="Volcano Plot")


  # top KM #
  top.genes <- unique(gsub(" ","",gsub("-1\\.5| 2| \\-2","",as.character(uni$Feature))))[1:10]
  top.genes <- top.genes[!is.na(top.genes)]
  if(any(apply(X %>% select(top.genes),2,
               function(x){length(unique(x))/length(x)}) > 0.1)){
    rm <- as.numeric(which(apply(X %>% select(top.genes),2,
                                 function(x){length(unique(x))/length(x)}) > 0.1))
    warning(paste0("Some covariate were seemingly continuous and therefore the Kaplan-Meier
            estimates will not be calculated for the following:",
                   paste0(colnames(X %>% select(top.genes))[rm],collapse = ",")))
    X <- (X %>% select(top.genes))[,-rm]
    top.genes <- top.genes[-rm]
  }

  KM.plots <- lapply(top.genes,function(x){
    y <- X[,x]
    if(length(unique(y[!is.na(y)]))==2)
      y <- factor(ifelse(X[,x] == 1,"Mutant","WildType"),levels = c("WildType","Mutant"))
    else
      y <- factor(as.numeric(as.character(y)),
                  levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(y)))])
    temp <- as.data.frame(cbind(surv.dat,y))
    if(LT == FALSE) fit <- survfit(Surv(time,status)~y,data=temp)
    if(LT == TRUE) fit <- survfit(Surv(time1,time2,status)~y,data=temp)

    if(length(unique(y[!is.na(y)]))==2)
      try(ggsurvplot(
        fit,
        data = temp,
        size = 1,
        palette =
          c("#E7B800", "#2E9FDF"),
        conf.int = TRUE,
        pval = ifelse(LT,FALSE,TRUE),
        risk.table = TRUE,
        legend.labs =
          c("WildType", "Mutant"),
        risk.table.col = "strata",
        risk.table.height = 0.25,
        ggtheme = theme_bw()
      ) + labs(title = x),silent = TRUE)

    else
      try(ggsurvplot(
        fit,
        data = temp,
        size = 1,
        conf.int = TRUE,
        pval = ifelse(LT,FALSE,TRUE),
        risk.table = TRUE,
        risk.table.col = "strata",
        risk.table.height = 0.25,
        ggtheme = theme_bw()
      ) + labs(title = x),silent = TRUE)
  })

  return(list("tab" = uni,"p"=uniVolcano,"KM"=KM.plots))
}

