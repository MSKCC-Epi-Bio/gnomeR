#' gen_summary
#'
#' Creates a summary table of the distribution of the genetic features by a specific outcome/covariate of interest. The difference
#' is tested using Fisher's exact test and further adjusted for multiple comparisons. Note that continuous genetic factors
#' are dichotomized at their median.
#'
#' @param gen_dat A matrix or dataframe, with samples as rows and features as columns.
#' @param outcome A leveled vector of length equal to the number of rows in gen_dat.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of samples
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @param cont Should the outcome be treated as a continuous value. Default is FALSE treated as categorical.
#' @param rank Should the table returned be ordered by Pvalue. Boolean, default is T
#' @param na_filter A numeric value between 0 and 1 (1 not included) that is the upper bound for the proportion of missing
#' values in the features of the inputted gen_dat matrix. Variables that exceed this proportion of missing values will be removed.
#' @return fits : a table of odds ratio and pvalues.
#' @return forest.plot : A forest plot of the top 10 hits.
#' @export
#' @examples library(gnomeR)
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))
#' ## binary outcome ##
#' outcome <- as.character(clin.sample$Sample.Type[match(samples,clin.sample$Sample.Identifier)])
#' gen_dat <- create_gene_binary(samples = samples,mutation = mut)

#' ## Continuous outcome ##
#' set.seed(1)
#' outcome <-  rnorm(n = nrow(gen_dat))
#' tab.out <- gen_summary(gen_dat = gen_dat,
#'                    outcome = outcome,
#'                    filter = 0.05,
#'                    cont = TRUE, rank = TRUE)
#' tab.out$fits
#' tab.out$vPlot
#' @import
#' dplyr
#' tibble


gen_summary <- function (gen_dat, outcome, filter = 0, cont = F, rank = T, na_filter = 0){

  # perform checks #
  if(filter < 0 || filter >= 1)
    stop("The filter should be between 0 and 1 (1 non included) to proceed.")

  if(na_filter < 0 || na_filter >= 1)
    stop("The filter for missing proportion should be between 0 and 1 (1 non included) to proceed.")

  if (length(unique(outcome))/length(outcome) > 0.5 && cont == F){
    warning("The outcome you provided had too many unique values, and will be considered continuous.")
    cont = T
    outcome <- as.numeric(as.character(outcome))
  }
  # check if there are any columns that are all NAs or have a single possible value #
  if (length(which(apply(gen_dat, 2, function(x) {
    length(unique(x[!is.na(x)])) <= 1
  } || all(is.na(x)) ))) > 0)
    # if found then remove #
    gen_dat <- gen_dat[, -which(apply(gen_dat, 2, function(x) {
      length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))
    }))]


  # apply filter to remove rare events #
  if (filter > 0) {

    # find those to remove #
    rm <- apply(gen_dat, 2, function(x) {
      any(summary(as.factor(x[!is.na(x)]))/length(x) < filter) # length(x[!is.na(x)])
    })
    genes.rm <- names(rm[which(rm)])
    gen_dat <- gen_dat %>% select(-one_of(genes.rm))
  }

  # apply filter for missing values #
  if(na_filter > 0){
    rm <- apply(gen_dat, 2, function(x) {
      sum(is.na(x))/length(x) > na_filter
    })
    genes.rm <- names(rm[which(rm)])
    if(length(genes.rm) > 0)
      gen_dat <- gen_dat %>% select(-one_of(genes.rm))
  }

  # check that some features are left after filtering #
  if (is.null(dim(gen_dat)) || dim(gen_dat)[2] == 0)
    stop("Only one or fewer genes are left after filtering. We need a minimum of two. Please relax the filter argument.")


  # Categorical outcomes tests through fisher test #
  if (!cont) {
    # make sure the outcome is a factor #
    if (!is.factor(outcome))
      outcome <- as.factor(outcome)
    # counter <<- 0

    fits <- lapply(gen_dat, function(x){
      # counter <<- counter + 1
      # print(counter)
      # print(x)
      # if a feature is continuous split it at the median #
      if (length(unique(x))/length(x) > 0.5) {
        x <- as.numeric(as.character(x))
        x <- ifelse(x > stats::median(x, na.rm = T), 1, 0)
      }

      # perform test #
      test <- stats::fisher.test(x, outcome)
      # print(is.numeric(x))
      # for binary data #
      if(is.numeric(x)){
        out <- c(paste0(round(sum(x[!is.na(x)])/length(x)*100,3) ,"%"))
        for(i in levels(outcome)){
          out <- c(out,paste0(round(sum(x[!is.na(x) & outcome == i])/sum(outcome == i)*100,3) ,"%"))
        }
        out <- c(out,
                 ifelse(!is.null(test$estimate),round(test$estimate,2), ""),
                 ifelse(!is.null(test$p.value),formatC(test$p.value, format = "e",digits = 2), ""),
                 ifelse(!is.null(test$conf.int[1]),round(test$conf.int[1], 2), ""),
                 ifelse(!is.null(test$conf.int[2]),round(test$conf.int[2], 2), "")
        )
        names(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
                        "OddsRatio","Pvalue", "Lower", "Upper")
      }

      # for categorical data #
      if(is.factor(x)){
        out <- as.data.frame(paste0(round(c(as.numeric(summary(x[!is.na(x)]))/length(x))*100,3) ,"%"))
        colnames(out) <- "Overall"
        rownames(out) <- levels(x)
        for(i in levels(outcome)){
          add.out <- c()
          for(j in levels(x)){
            add.out <- c(add.out,paste0(round(sum(x[outcome == i] == j,na.rm = T)/sum(outcome == i)*100,3),"%"))
          }
          out[,i] <- add.out
        }
        out$OddsRatio <- c("",ifelse(!is.null(test$estimate),round(test$estimate,2), ""),
                           rep("",nrow(out)-2))
        out$Pvalue <- c("",ifelse(!is.null(test$p.value),formatC(test$p.value, format = "e",digits = 2), ""),
                        rep("",nrow(out)-2))
        out$Lower <- c("",ifelse(!is.null(test$conf.int[1]),round(test$conf.int[1], 2), ""),
                       rep("",nrow(out)-2))
        out$Upper <- c("",ifelse(!is.null(test$conf.int[2]),round(test$conf.int[2], 2), ""),
                       rep("",nrow(out)-2))
        out <- out[-1,]

      }

      return(out)
    })

    ps <- c()
    for(i in 1:length(fits)){
      if(is.null(dim(fits[[i]])))
        ps[i] <- fits[[i]]["Pvalue"]
      else{
        ps[i] <- fits[[i]][1,"Pvalue"]
        if(nrow(fits[[i]])==1)
          names(fits)[[i]] <- paste0(names(fits)[[i]],".",rownames(fits[[i]]))
      }
    }
    if(rank)
      fits <- as.data.frame(do.call('rbind',fits[order(as.numeric(ps))]))
    else
      fits <- as.data.frame(do.call('rbind',fits))

    fits$Feature <- rownames(fits)
    fits$FDR <- stats::p.adjust(p = as.numeric(fits$Pvalue), method = "fdr")
    fits <- fits[,c("Feature","Overall", levels(outcome), "OddsRatio", "Lower",
                    "Upper", "Pvalue", "FDR")]
    colnames(fits)[3:(length(levels(outcome)) + 2)] <- paste0(colnames(fits)[3:(length(levels(outcome)) +
                                                                                  2)], "(N=", as.numeric(summary(outcome)), ")")

    if (any(as.character(fits$OddsRatio) != "")) {
      f.dat <- fits
      f.dat <- f.dat %>%
        mutate(Pvalue = as.numeric(as.character(.data$Pvalue)),
               OddsRatio = as.numeric(as.character(.data$OddsRatio))) %>%
        filter(!is.infinite(.data$OddsRatio)) %>%
        mutate(OddsRatio = as.numeric(as.character(.data$OddsRatio)),
               Lower = as.numeric(as.character(.data$Lower)), Upper = as.numeric(as.character(.data$Upper)),
        )
      f.dat <- f.dat[1:min(10, nrow(f.dat)), ]

      forest.plot <- f.dat %>%
        ggplot(aes(x = .data$Feature, y = .data$OddsRatio,
                   ymin = .data$Lower, ymax = .data$Upper)) + geom_pointrange(aes(col = .data$Feature)) +
        geom_hline(aes(fill = .data$Feature), yintercept = 1,
                   linetype = 2) + xlab("Feature") + ylab("Risk Ratio (95% Confidence Interval)") +
        geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper,
                          col = .data$Feature), width = 0.5, cex = 1) + coord_flip()

      vplot <- plot_ly(data = fits %>%
                         tibble::rownames_to_column("GeneName") %>%
                         mutate(Pvalue = as.numeric(as.character(.data$Pvalue)),
                                OddsRatio = as.numeric(as.character(.data$OddsRatio)),
                                FDRsign = ifelse(as.numeric(as.character(.data$FDR)) <
                                                   0.05, "Significant", "Non signifcant")) %>%
                         filter(!is.na(.data$Pvalue),
                                is.numeric(.data$Pvalue)), x = ~.data$OddsRatio,
                       y = ~-log10(.data$Pvalue),color = .data$PFDRsign,
                       text = ~paste("Gene :", GeneName,
                                     "</br> Odds Ratio :",
                                     round(.data$OddsRatio, digits = 2)), mode = "markers") %>%
        layout(title = "Volcano Plot")
    }
    else {
      forest.plot <- NULL
      vplot <- NULL
    }

    return(list(fits = fits, forest.plot = forest.plot, vPlot = vplot))
  }


  if (cont) {

    fits <- lapply(gen_dat,function(x){
      fit <- stats::lm(outcome ~ x)
      temp <- as.data.frame(summary(fit)$coefficient)
      colnames(temp) <- c("Estimate", "SD", "tvalue",
                          "Pvalue")
      temp$Lower <- round(temp$Estimate - 1.96*temp$SD,3)
      temp$Upper <- round(temp$Estimate + 1.96*temp$SD,3)
      if(is.numeric(x))
        temp$EventFrequency <- round(sum(x[!is.na(x)])/length(x),digits = 3)
      else
        temp$EventFrequency <- round(summary(x[!is.na(x)])/length(x),digits = 3)
      temp <- temp[-1,]
      return(temp)
    })

    if(rank){
      ps <- c()
      for(i in 1:length(fits)){
        ps[i] <- min(fits[[i]]$Pvalue)
        if(nrow(fits[[i]])==1)
          names(fits)[[i]] <- paste0(names(fits)[[i]],".",rownames(fits[[i]]))
      }
      fits <- as.data.frame(do.call('rbind',fits[order(as.numeric(ps))]))
    }
    else
      fits <- as.data.frame(do.call('rbind',fits))

    fits <- fits %>%
      mutate(Feature = gsub("\\.","_",gsub("\\.x"," ",rownames(fits))),
             FDR = formatC(stats::p.adjust(.data$Pvalue, method = 'fdr'), format="e", digits=2),
             Pvalue = formatC(.data$Pvalue, format="e", digits=2),
             Estimate = round(.data$Estimate, 3),
             SD = round(.data$SD, 3)
      ) %>%
      select(.data$Feature, .data$EventFrequency, .data$Estimate,
             .data$SD, .data$Lower, .data$Upper, .data$Pvalue, .data$FDR)


    vPlot <- plot_ly(data = fits %>%
                       mutate(FDRsign = ifelse(as.numeric(as.character(.data$FDR)) <
                                                 0.05, "Significant", "Non signifcant"),
                              Pvalue = as.numeric(.data$Pvalue)) %>%
                       filter(!is.na(.data$Pvalue),
                              is.numeric(.data$Pvalue)), x = ~.data$Estimate,
                     y = ~-log10(.data$Pvalue),color = ~FDRsign,
                     text = ~paste("Feature :", Feature, "</br> Estimate :",
                                   round(.data$Estimate, digits = 2)), mode = "markers",
                     size = ~.data$EventFrequency) %>%
      layout(title = "Volcano Plot",
             xaxis = list(title = "Estimate"),
             yaxis = list(title = "-log10(Pvalue)"))


    return(list(fits = fits, vPlot = vPlot))
  }
}


