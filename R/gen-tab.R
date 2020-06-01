#' gen.tab
#'
#' Creates a summary table of the distribution of the genetic features by a specific outcome/covariate of interest. The difference
#' is tested using Fisher's exact test and further adjusted for multiple comparisons. Note that continuous genetic factors
#' are dichotomized at their median.
#'
#' @param gen.dat A matrix or dataframe, with patients as rows and features as columns.
#' @param outcome A leveled vector of length equal to the number of rows in gen.dat.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of patients
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @param paired Boolean if the data are paired. Default is FALSE.
#' @param cont Should the outcome be treated as a continuous value. Default is FALSE treated as categorical.
#' @param rank Should the table returned be ordered by Pvalue. Boolean, default is T
#' @return fits : a table of odds ratio and pvalues.
#' @return forest.plot : A forest plot of the top 10 hits.
#' @export
#' @examples library(gnomeR)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))
#' ## binary outcome ##
#' outcome <- as.character(clin.sample$Sample.Type[match(patients,clin.sample$Sample.Identifier)])
#' gen.dat <- binmat(patients = patients,maf = mut)
#' gen.tab(gen.dat = gen.dat,
#'         outcome = outcome,
#'         filter = 0.05,paired = FALSE,
#'         cont = FALSE,rank = TRUE)
#' ## Continuous outcome ##
#' set.seed(1)
#' outcome <-  rnorm(n = nrow(gen.dat))
#' tab.out <- gen.tab(gen.dat = gen.dat,
#'                    outcome = outcome,
#'                    filter = 0.05,paired = FALSE,
#'                    cont = TRUE,rank = TRUE)
#' tab.out$fits
#' tab.out$vPlot
#' @import
#' ggrepel
#' exact2x2

gen.tab <- function (gen.dat, outcome, filter = 0, paired = F, cont = F,
                     rank = T)
{
  if(filter < 0 || filter >= 1)
    stop("The filter should be between 0 and 1 (1 non included) to proceed.")

  if (length(which(apply(gen.dat, 2, function(x) {
    length(unique(x[!is.na(x)])) <= 1
  } || all(is.na(x))))) > 0)
    gen.dat <- gen.dat[, -which(apply(gen.dat, 2, function(x) {
      length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))
    }))]
  if (filter > 0) {

    rm <- apply(gen.dat, 2, function(x) {
      # if(is.numeric(x))
      #   sum(x, na.rm = T)/length(x) < filter
      # else
      any(summary(as.factor(x[!is.na(x)]))/length(x[!is.na(x)]) < filter)
    })
    genes.rm <- names(rm[which(rm)])
    length(genes.rm)
    gen.dat <- gen.dat %>% select(-one_of(genes.rm))
  }
  if (is.null(dim(gen.dat)) || dim(gen.dat)[2] == 0)
    stop("Only one or fewer genes are left after filtering. We need a minimum of two. Please relax the filter argument.")
  if (!cont) {
    if (!is.factor(outcome))
      outcome <- as.factor(outcome)

    fits <- lapply(gen.dat, function(x) {
      if (length(unique(x))/length(x) > 0.5) {
        x <- as.numeric(as.character(x))
        x <- ifelse(x > stats::median(x, na.rm = T), 1, 0)
      }
      else if(length(unique(x[!is.na(x)])) > 2)
        x <- factor(x,
                    levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(x)))])
      else if(length(unique(x[!is.na(x)])) == 2)
        x <- as.numeric(as.character(x))

      if (paired == F)
        test <- stats::fisher.test(x, outcome)
      if (paired == T) {
        test <- mcnemar.exact(x = x[1:(length(x)/2)],
                              y = x[(length(x)/2 + 1):length(x)])
      }
      if(length(unique(x[!is.na(x)])) == 2)
        out <- c(sum(x, na.rm = T)/length(x))
      else
        out <- summary(x)/length(x)
      for (i in 1:length(levels(outcome))) {
        if(length(unique(x[!is.na(x)])) == 2)
          out <- c(out, sum(x[which(outcome == levels(outcome)[i])],
                            na.rm = T)/length(which(outcome == levels(outcome)[i])))
        else
          out <- cbind(out,summary(x[which(outcome == levels(outcome)[i])])/length(which(outcome == levels(outcome)[i])))
      }
      if(is.null(dim(out)))
        out <- paste0(round(as.numeric(out) * 100, digits = 2),
                      "%")
      else
        out <- apply(out,2, function(x){
          paste0(round(as.numeric(x) * 100, digits = 2),
                 "%")
        })
      if (!is.null(test$estimate)) {
        out <- c(out, round(as.numeric(test$estimate),
                            digits = 2),
                 formatC(test$p.value, format = "e",
                         digits = 2),
                 formatC(stats::p.adjust(test$p.value, method ="fdr",n = ncol(gen.dat)), format = "e",
                         digits = 2),
                 round(as.numeric(test$conf.int),
                       digits = 2))
        names(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
                        "OddsRatio", "Pvalue", "FDR", "Lower", "Upper")
      }
      else {
        if(length(unique(x[!is.na(x)])) == 2){
          out <- c(out, "",
                   formatC(test$p.value, format = "e",
                           digits = 2),
                   formatC(stats::p.adjust(test$p.value, method ="fdr",n = ncol(gen.dat)), format = "e",
                           digits = 2),
                   "","")
          names(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
                          "OddsRatio", "Pvalue", "FDR", "Lower", "Upper")
        }
        else{
          out <- cbind(out, rep("",nrow(out)),
                       c(formatC(test$p.value, format = "e",
                                 digits = 2),rep("",nrow(out)-1)),
                       formatC(stats::p.adjust(test$p.value, method ="fdr",n = ncol(gen.dat)), format = "e",
                               digits = 2),
                       rep("",nrow(out)),rep("",nrow(out)))
          colnames(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
                             "OddsRatio","Pvalue", "FDR", "Lower", "Upper")
        }
      }
      return(out)
    })

    if(rank){
      ps <- c()
      for(i in 1:length(fits)){
        if(is.null(dim(fits[[i]])))
          ps[i] <- fits[[i]][5]
        else{
          ps[i] <- fits[[i]][1,5]
          rownames(fits[[i]]) <- c(names(fits)[[i]],rep("",nrow(fits[[i]])-1))
        }
      }
      fits <- do.call('rbind',fits[order(as.numeric(ps))])
    }
    else{
      for(i in 1:length(temp)){
        if(!is.null(dim(temp[[i]])))
          rownames(temp[[i]]) <- c(names(temp)[[i]],rep("",nrow(temp[[i]])-1))
      }
      fits <- do.call('rbind',fits)
    }

    fits <- as.data.frame(cbind(rownames(fits),fits))
    colnames(fits)[1] <- "Feature"
    colnames(fits)[3:(length(levels(outcome)) + 2)] <- paste0(colnames(fits)[3:(length(levels(outcome)) +
                                                                                  2)], "(N=", as.numeric(summary(outcome)), ")")
    fits$FDR <- formatC(stats::p.adjust(as.numeric(as.character(fits$Pvalue)),
                                        method = "fdr"), format = "e", digits = 2)
    if (rank)
      fits <- fits[order(as.numeric(as.character(fits$Pvalue))),
                   ]
    if (any(as.character(fits$OddsRatio) != "")) {
      f.dat <- fits
      f.dat <- f.dat %>% filter(!is.infinite(.data$OddsRatio)) %>%
        mutate(OddsRatio = as.numeric(as.character(.data$OddsRatio)),
               Lower = as.numeric(as.character(.data$Lower)), Upper = as.numeric(as.character(.data$Upper)),
        )
      f.dat <- f.dat[1:min(10, nrow(f.dat)), ]
      forest.plot <- f.dat %>% ggplot(aes(x = .data$Feature, y = .data$OddsRatio,
                                          ymin = .data$Lower, ymax = .data$Upper)) + geom_pointrange(aes(col = .data$Feature)) +
        geom_hline(aes(fill = .data$Feature), yintercept = 1,
                   linetype = 2) + xlab("Feature") + ylab("Risk Ratio (95% Confidence Interval)") +
        geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper,
                          col = .data$Feature), width = 0.5, cex = 1) + coord_flip()
      vplot <- fits  %>%
        mutate(Pvalue = as.numeric(as.character(.data$Pvalue)),
               OddsRatio = as.numeric(as.character(.data$OddsRatio)),
               FDRsign = ifelse(as.numeric(as.character(.data$FDR)) <
                                  0.05, "Significant", "Non signifcant")) %>%
        filter(!is.na(.data$Pvalue)) %>%
        ggplot(aes(x = .data$OddsRatio, y = -log10(.data$Pvalue),
                   fill = .data$FDRsign, color = .data$FDRsign)) + geom_point() +
        geom_label_repel(aes(label = ifelse(.data$FDRsign ==
                                              "Significant", as.character(.data$Feature), "")), color = "white")
    }
    else {
      forest.plot <- NULL
      vplot <- NULL
    }

    return(list(fits = fits, forest.plot = forest.plot, vPlot = vplot))
  }
  if (cont) {
    fits <- as.data.frame(do.call("rbind", apply(gen.dat,
                                                 2, function(x) {
                                                   if (length(unique(x[!is.na(x)]))/length(x[!is.na(x)]) > 0.25 ||
                                                       length(unique(x[!is.na(x)])) == 2) {
                                                     x <- as.numeric(x)
                                                     fit <- stats::lm(outcome ~ x)
                                                     temp <- as.data.frame(summary(fit)$coefficient)
                                                     colnames(temp) <- c("Estimate", "SD", "tvalue",
                                                                         "Pvalue")
                                                     if (is.numeric(x))
                                                       temp$MutationFreq <- sum(x, na.rm = T)/length(x[!is.na(x)])
                                                     else temp$MutationFreq <- 0
                                                     out <- temp[2, c(1, 2, 4, 5)]
                                                   }
                                                   else {
                                                     x <- factor(x,
                                                                 levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in%
                                                                                                         as.numeric(as.character(x)))])
                                                     fit <- stats::lm(outcome ~ x)
                                                     temp <- as.data.frame(summary(fit)$coefficient)
                                                     colnames(temp) <- c("Estimate", "SD", "tvalue",
                                                                         "Pvalue")
                                                     temp$MutationFreq <- summary(x)/length(x)
                                                     out <- as.data.frame(temp[2:nrow(temp), c(1,
                                                                                               2, 4, 5)])
                                                     rownames(out) <- gsub("x", "",
                                                                           rownames(out))
                                                     # rownames(out) <- gsub("as.factor\\(x\\)", "",
                                                     #                       rownames(out))
                                                   }
                                                   return(out)
                                                 })))
    fits$FDR <- stats::p.adjust(fits$Pvalue, method = "fdr")
    fits$GeneName <- rownames(fits)
    if (all(apply(gen.dat, 2, is.numeric))) {

      vPlot <- try(plot_ly(data = fits %>% filter(!is.na(.data$Pvalue),
                                                  is.numeric(.data$Pvalue)), x = ~Estimate, y = ~-log10(.data$Pvalue),
                           text = ~paste("Gene :", GeneName, "</br> Estimate :",
                                         round(Estimate, digits = 2)), mode = "markers",
                           size = ~MutationFreq, color = ~Estimate) %>%
                     layout(title = "Volcano Plot"), silent = T)
    }
    else {
      vPlot <- try(plot_ly(data = fits %>% filter(!is.na(.data$Pvalue),
                                                  is.numeric(.data$Pvalue)), x = ~Estimate, y = ~-log10(.data$Pvalue),
                           text = ~paste("Gene :", GeneName, "</br> Estimate :",
                                         round(Estimate, digits = 2)), mode = "markers") %>%
                     layout(title = "Volcano Plot"), silent = T)
    }
    if (rank)
      fits <- fits[order(fits$Pvalue), -ncol(fits)]
    return(list(fits = fits, vPlot = vPlot))
  }
}


# function (gen.dat, outcome, filter = 0, paired = F, cont = F,
#                      rank = T)
# {
#   if(filter < 0 || filter >= 1)
#     stop("The filter should be between 0 and 1 (1 non included) to proceed.")
#
#   if (length(which(apply(gen.dat, 2, function(x) {
#     length(unique(x[!is.na(x)])) == 1
#   } || all(is.na(x))))) > 0)
#     gen.dat <- gen.dat[, -which(apply(gen.dat, 2, function(x) {
#       length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))
#     }))]
#   if (filter > 0) {
#     temp <- apply(gen.dat, 2, function(x) {
#       length(unique(x[!is.na(x)])) == 2
#     })
#     genes.bin <- names(temp[which(temp)])
#     if (length(genes.bin) == ncol(gen.dat))
#       rm <- apply(gen.dat, 2, function(x) {
#         if(is.numeric(x))
#           sum(x, na.rm = T)/length(x) < filter
#         else
#           all(summary(as.factor(x[!is.na(x)]))/length(x) > filter)
#       })
#     else rm <- apply(gen.dat[, genes.bin], 2, function(x) {
#       if(is.numeric(x))
#         sum(x, na.rm = T)/length(x) < filter
#       else
#         all(summary(as.factor(x[!is.na(x)]))/length(x) > filter)
#     })
#     genes.rm <- names(rm[which(rm)])
#     gen.dat <- gen.dat %>% select(-one_of(genes.rm))
#   }
#   if (is.null(dim(gen.dat)) || dim(gen.dat)[2] == 0)
#     stop("Only one or fewer genes are left after filtering. We need a minimum of two. Please relax the filter argument.")
#   if (!cont) {
#     if (is.character(outcome))
#       outcome <- as.factor(outcome)
#     fits <- as.data.frame(t(apply(gen.dat, 2, function(x) {
#       if (length(unique(x))/length(x) > 0.5) {
#         x <- ifelse(x > stats::median(x, na.rm = T), 1, 0)
#       }
#       if (paired == F)
#         test <- stats::fisher.test(x, outcome)
#       if (paired == T) {
#         test <- mcnemar.exact(x = x[1:(length(x)/2)],
#                               y = x[(length(x)/2 + 1):length(x)])
#       }
#       out <- c(sum(x, na.rm = T)/length(x))
#       for (i in 1:length(levels(outcome))) {
#         out <- c(out, sum(x[which(outcome == levels(outcome)[i])],
#                           na.rm = T)/length(which(outcome == levels(outcome)[i])))
#       }
#       out <- paste0(round(as.numeric(out) * 100, digits = 2),
#                     "%")
#       if (!is.null(test$estimate)) {
#         out <- c(out, round(as.numeric(test$estimate),
#                             digits = 2), formatC(test$p.value, format = "e",
#                                                  digits = 2), round(as.numeric(test$conf.int),
#                                                                     digits = 2))
#         names(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
#                         "OddsRatio", "Pvalue", "Lower", "Upper")
#       }
#       else {
#         out <- c(out, test$p.value, round(as.numeric(test$conf.int),
#                                           digits = 2))
#         names(out) <- c("Overall", levels(outcome)[1:length(levels(outcome))],
#                         "Pvalue")
#       }
#       return(out)
#     })))
#     colnames(fits)[2:(length(levels(outcome)) + 1)] <- paste0(colnames(fits)[2:(length(levels(outcome)) +
#                                                                                   1)], "(N=", as.numeric(summary(outcome)), ")")
#     fits$FDR <- formatC(stats::p.adjust(as.numeric(as.character(fits$Pvalue)),
#                                         method = "fdr"), format = "e", digits = 2)
#     if (rank)
#       fits <- fits[order(as.numeric(as.character(fits$Pvalue))),
#                    ]
#     if (!is.null(fits$OddsRatio)) {
#       f.dat <- fits
#       f.dat$Gene <- rownames(f.dat)
#       f.dat <- f.dat %>% filter(!is.infinite(.data$OddsRatio)) %>%
#         mutate(OddsRatio = as.numeric(as.character(.data$OddsRatio)),
#                Lower = as.numeric(as.character(.data$Lower)), Upper = as.numeric(as.character(.data$Upper)),
#         )
#       f.dat <- f.dat[1:min(10, nrow(f.dat)), ]
#       forest.plot <- f.dat %>% ggplot(aes(x = .data$Gene, y = .data$OddsRatio,
#                                           ymin = .data$Lower, ymax = .data$Upper)) + geom_pointrange(aes(col = .data$Gene)) +
#         geom_hline(aes(fill = .data$Gene), yintercept = 1,
#                    linetype = 2) + xlab("Gene") + ylab("Risk Ratio (95% Confidence Interval)") +
#         geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper,
#                           col = .data$Gene), width = 0.5, cex = 1) + coord_flip()
#       vplot <- fits %>% rownames_to_column("Gene") %>%
#         mutate(Pvalue = as.numeric(as.character(.data$Pvalue)),
#                OddsRatio = as.numeric(as.character(.data$OddsRatio)),
#                FDRsign = ifelse(as.numeric(as.character(.data$FDR)) <
#                                   0.05, "Significant", "Non signifcant")) %>%
#         ggplot(aes(x = .data$OddsRatio, y = -log10(.data$Pvalue),
#                    fill = .data$FDRsign, color = .data$FDRsign)) + geom_point() +
#         geom_label_repel(aes(label = ifelse(.data$FDRsign ==
#                                               "Significant", as.character(.data$Gene), "")), color = "white")
#     }
#     else {
#       forest.plot <- NULL
#       vplot <- NULL
#     }
#     return(list(fits = fits, forest.plot = forest.plot, vPlot = vplot))
#   }
#   if (cont) {
#     fits <- as.data.frame(do.call("rbind", apply(gen.dat,
#                                                  2, function(x) {
#                                                    if (length(unique(x))/length(x) > 0.25 || length(unique(x)) ==
#                                                        2) {
#                                                      fit <- stats::lm(outcome ~ x)
#                                                      temp <- as.data.frame(summary(fit)$coefficient)
#                                                      colnames(temp) <- c("Estimate", "SD", "tvalue",
#                                                                          "Pvalue")
#                                                      if (is.numeric(x))
#                                                        temp$MutationFreq <- sum(x, na.rm = T)/length(x[!is.na(x)])
#                                                      else temp$MutationFreq <- 0
#                                                      out <- temp[2, c(1, 2, 4, 5)]
#                                                    }
#                                                    else {
#                                                      fit <- stats::lm(outcome ~ as.factor(x))
#                                                      temp <- as.data.frame(summary(fit)$coefficient)
#                                                      colnames(temp) <- c("Estimate", "SD", "tvalue",
#                                                                          "Pvalue")
#                                                      temp$MutationFreq <- 0
#                                                      out <- as.data.frame(temp[2:nrow(temp), c(1,
#                                                                                                2, 4, 5)])
#                                                      rownames(out) <- gsub("as.factor\\(x\\)", "",
#                                                                            rownames(out))
#                                                    }
#                                                    return(out)
#                                                  })))
#     fits$FDR <- stats::p.adjust(fits$Pvalue, method = "fdr")
#     fits$GeneName <- rownames(fits)
#     if (all(apply(gen.dat, 2, is.numeric))) {
#
#       vPlot <- try(plot_ly(data = fits %>% filter(!is.na(.data$Pvalue),
#                                                   is.numeric(.data$Pvalue)), x = ~Estimate, y = ~-log10(.data$Pvalue),
#                            text = ~paste("Gene :", GeneName, "</br> Estimate :",
#                                          round(Estimate, digits = 2)), mode = "markers",
#                            size = ~MutationFreq, color = ~Estimate) %>%
#                      layout(title = "Volcano Plot"), silent = T)
#     }
#     else {
#       vPlot <- try(plot_ly(data = fits %>% filter(!is.na(.data$Pvalue),
#                                                   is.numeric(.data$Pvalue)), x = ~Estimate, y = ~-log10(.data$Pvalue),
#                            text = ~paste("Gene :", GeneName, "</br> Estimate :",
#                                          round(Estimate, digits = 2)), mode = "markers") %>%
#                      layout(title = "Volcano Plot"), silent = T)
#     }
#     if (rank)
#       fits <- fits[order(fits$Pvalue), -ncol(fits)]
#     return(list(fits = fits, vPlot = vPlot))
#   }
# }
