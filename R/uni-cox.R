#' uni.cox
#' Performs univariate cox proportional hazard model on every feature
#' @param X Matrix/surv.datframe of genomic features, continuous or binary (note cannot handle categorical surv.dat for the moment).
#' @param surv.dat a surv.dat frame containing the survival information. This can be made of 2 or 3 columns. 1 or 2 for time,
#' and one for status (where 1 is event and 0 is no event).
#' @param surv.formula a survival formula with names matching those in surv.dat eg: Surv(time,status)~.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of samples
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @param genes a character vector of gene names that will be the only ones to be kept. Default is NULL, all genes are used.
#' @param na.filt A numeric value between 0 and 1 (1 not included) that is the upper bound for the proportion of missing
#' values in the features of the inputted gen.dat matrix. Variables that exceed this proportion of missing values will be removed.
#' @return tab A table of all the fits performed sorted by adjusted pvalues.
#' @return p An interactive plot of log(pvalue) by hazard ration.
#' @return KM List of survival plots of the top 10 most significant genes
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' gen.dat <- binary_matrix(samples = samples, mutation = mut)
#' surv.dat <- clin.patients %>%
#' filter(X.Patient.Identifier %in%
#' abbreviate(samples,strict = TRUE, minlength = 9)) %>%
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


uni.cox <- function(X,surv.dat,surv.formula,filter = 0,genes = NULL, na.filt = 0){

  # filtering #
  if(!(filter >= 0 && filter < 1))
    stop("Please select a filter value between 0 and 1")
  if(na.filt < 0 || na.filt >= 1)
    stop("The filter for missing proportion should be between 0 and 1 (1 non included) to proceed.")

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
  if(length(which(apply(X, 2, function(x){length(unique(x[!is.na(x)])) == 1} || all(is.na(x))))) > 0){
    X <- X[,-which(apply(X, 2, function(x){length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))}))]
    # relevel cna data #
    if(length(grep("\\.cna",colnames(X))) > 0){
      X <- cbind(X %>%
                   select(-c(ends_with(".cna"))),
                 X %>%
                   select(ends_with(".cna")) %>%
                   mutate_all(
                     ~factor(.,
                             levels = c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION")[
                               which(c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION") %in% .)
                               ])
                   )
      )
    }
  }

  # apply filter to remove rare events #
  if (filter > 0) {

    # find those to remove #
    rm <- apply(X, 2, function(x) {
      any(summary(as.factor(x[!is.na(x)]))/length(x) < filter) # length(x[!is.na(x)])
    })
    genes.rm <- names(rm[which(rm)])
    X <- X %>% select(-one_of(genes.rm))
  }

  # apply filter for missing values #
  if(na.filt > 0){
    rm <- apply(X, 2, function(x) {
      sum(is.na(x))/length(x) > na.filt
    })
    genes.rm <- names(rm[which(rm)])
    if(length(genes.rm) > 0)
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
    datSurv <- with(surv.dat,Surv(time,status))
  }
  if(length(as.list(survResponse)) == 4){
    colnames(surv.dat)[match(as.list(survResponse)[2:4],colnames(surv.dat))] <- c("time1","time2","status")
    LT = TRUE
    timevars <- match(c("time1","time2","status"),colnames(surv.dat))
    surv.dat <- surv.dat[,c(timevars,
                            c(1:ncol(surv.dat))[-timevars])]
    datSurv <- with(surv.dat,Surv(time1,time2,status))
  }


  fits <- lapply(X, function(x){
    fit <- coxph(datSurv ~ x)
    temp <- as.data.frame(summary(fit)$coefficient)
    colnames(temp) <- c("Estimate", "HR", "SD", "tvalue",
                        "Pvalue")
    temp$Lower <- round(exp(temp$Estimate - 1.96*temp$SD),3)
    temp$Upper <- round(exp(temp$Estimate + 1.96*temp$SD),3)
    if(is.numeric(x))
      temp$EventFrequency <- round(sum(x[!is.na(x)])/length(x),digits = 3)
    else
      temp$EventFrequency <- round(summary(x[!is.na(x)])/length(x),digits = 3)[-1]
    return(temp)
  })

  ps <- c()
  for(i in 1:length(fits)){
    ps[i] <- min(fits[[i]]$Pvalue)
    if(nrow(fits[[i]])==1)
      names(fits)[[i]] <- paste0(names(fits)[[i]],".",rownames(fits[[i]]))
  }
  fits <- as.data.frame(do.call('rbind',fits[order(as.numeric(ps))]))

  fits <- fits %>%
    mutate(Feature = gsub("\\.","_",gsub("\\.x"," ",rownames(fits))),
           FDR = formatC(stats::p.adjust(.data$Pvalue, method = 'fdr'), format="e", digits=2),
           Pvalue = formatC(.data$Pvalue, format="e", digits=2),
           Estimate = round(.data$Estimate, 3),
           HR = round(.data$HR, 3),
           SD = round(.data$SD, 3)
    ) %>%
    select(.data$Feature, .data$EventFrequency, .data$Estimate,
           .data$HR, .data$SD, .data$Lower, .data$Upper,
           .data$Pvalue, .data$FDR)

  ### Volcano plot ###
  uniVolcano <- plot_ly(data = fits %>%
                          mutate(FDRsign = ifelse(as.numeric(as.character(.data$FDR)) <
                                                    0.05, "Significant", "Non signifcant"),
                                 Pvalue = as.numeric(.data$Pvalue)), x = ~.data$Estimate, y = ~-log10(.data$Pvalue),
                        text = ~paste('Feature :',.data$Feature,
                                      '<br> Hazard Ratio :',.data$HR,
                                      '<br> Event Frequency :',.data$EventFrequency),
                        mode = "markers",color = ~ifelse(.data$FDRsign == "Significant","blue","red")) %>%
    layout(title ="Volcano Plot")

  # top KM #
  top.genes <- gsub("_",".",unique(gsub(" .*","",as.character(fits$Feature)))[1:10])
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
    # if(length(unique(y[!is.na(y)]))==2)
    #   y <- factor(ifelse(X[,x] == 1,"Mutant","WildType"),levels = c("WildType","Mutant"))
    # else
    #   y <- factor(as.numeric(as.character(y)),
    #               levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(y)))])

    temp <- as.data.frame(cbind(surv.dat,y))
    if(LT == FALSE) fit <- survfit(Surv(time,status)~y,data=temp)
    if(LT == TRUE) fit <- survfit(Surv(time1,time2,status)~y,data=temp)

    ggsurvplot(
      fit,
      data = temp,
      size = 1,
      conf.int = TRUE,
      pval = ifelse(LT,FALSE,TRUE),
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.height = 0.25,
      ggtheme = theme_bw()
    ) + labs(title = x)

  })

  return(list("tab" = fits,"p"=uniVolcano,"KM"=KM.plots))
}

