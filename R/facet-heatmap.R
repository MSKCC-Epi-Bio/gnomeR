#' facets.heatmap
#'
#' Creates a heatmap of copy number profiles from segment files.
#'
#' @param seg a segmentation file containing the segmentation information of multiple samples
#' @param filenames the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param path the relative path to the files folder from your current directory
#' @param samples the names of the samples of the respective filenames. Default simply 1 to number of files.
#' @param min_purity the minimum purity of the sample required to be kept in the final dataset. Default is 0.3.
#' @param epsilon level of unions when aggregating segments between
#' @param ordered order in which samples should be printed. Default NUll leads to hierarchical clustering.
#' @param outcome for seg file only, if outcome associated with study it will be printed along the x axis for each patient
#' @param adaptive CNregions option to create adaptive segments
#' @return p a heatmap corresponding to the segment files inputted
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:1000]
#' samples.seg <- as.character(unlist(clin.sample %>%
#' filter(Sample.Identifier %in% samples,
#'  as.numeric(as.character(Tumor.Purity)) > 30) %>%
#'  select(Sample.Identifier)))
#' facet <- facets.heatmap(seg = seg, samples=samples.seg[0:100])
#' facet$p
#' @import
#' gplots
#' lattice
#' tibble


facets.heatmap <- function (seg = NULL, filenames = NULL, path = NULL, samples = NULL,
                            min_purity = 0.3, epsilon = 0.005, ordered = NULL, outcome = NULL,
                            adaptive = FALSE)
{
  if (is.null(seg) && is.null(filenames))
    stop("You must provide either a complete segmentation file\n         or a list of files to be loaded with their corresponding path")
  if (!is.null(seg) && !is.null(filenames))
    stop("Please provide either a complete segmentation file or a\n         list of segmentation files to be loaded")
  if (!is.null(filenames)) {
    dat <- facets_dat(seg = NULL, filenames, path, samples,
                      min_purity, epsilon, adaptive)
    reducedM <- dat$out.cn
    ploidy <- dat$ploidy
    purity <- dat$purity
    rownames(reducedM) <- as.character(abbreviate(rownames(reducedM),
                                                  minlength = 17))
    if(!is.null(outcome)) {
      missing <- which(is.na(match(names(outcome),rownames(reducedM))))
      if(length(missing)>0)
        outcome <- outcome[-missing]

      reducedM <- reducedM[match(names(outcome),rownames(reducedM)),]
    }

    imagedata = reducedM
    imagedata[imagedata > 1.5] = 1.5
    imagedata[imagedata < -1.5] = -1.5

    if (is.null(ordered)) {
      cl = stats::hclust(stats::dist(imagedata), method = "ward")
      imagedata.ordered = imagedata[cl$order, ]
      imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
    if (!is.null(ordered)) {
      if(length(missing) > 0)
        ordered <- order(outcome)
      imagedata.ordered = imagedata[ordered, ]
      imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
    chr = strsplit(colnames(imagedata), "\\.")
    chr = unlist(lapply(1:length(chr), function(x) chr[[x]][1]))
    chr = gsub("chr", "", chr)
    chr = as.numeric(chr)
    len = length(chr)
    chrom.ends <- rep(NA, length(table(chr)))
    d = 1
    for (r in unique(chr)) {
      chrom.ends[d] <- max(which(chr == r))
      d = d + 1
    }
    chrom.starts <- c(1, chrom.ends[-length(table(chr))] +
                        1)
    chrom.mids <- (chrom.starts + chrom.ends)/2
    bw = colorpanel(2, low = "white", high = "cadetblue4")
    colorkey = list(space = "right", height = 0.3, tick.number = 5)
    n <- nrow(reducedM)
    if (!is.null(outcome))
      x.lab <- outcome
    if (is.null(outcome))
      x.lab <- rep(" ", n)
    if (is.null(ordered))
      x.lab <- as.character(x.lab[cl$order])
    if (!is.null(ordered))
      x.lab <- as.character(x.lab[ordered])
    if (is.null(outcome) && is.null(ordered))
      scales = list(x = list(at = 1:n, labels = x.lab,
                             rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                    z = list(at = n:1, labels = purity[cl$order],
                             rot = 90))
    else scales = list(x = list(at = 1:n, labels = x.lab,
                                rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                       z = list(at = n:1, labels = rep(1, n), rot = 90))
    my.panel.levelplot.2 <- function(...) {
      panel.levelplot(...)
      panel.abline(h = len - chrom.starts[-1], col = "gray",
                   lwd = 1)
      panel.scales = list(x = list(at = 1:n), y = list(at = len -
                                                         chrom.mids), z = list())
    }
    my.panel = my.panel.levelplot.2
    p = levelplot(imagedata.ordered, panel = my.panel, scales = scales,
                  aspect = "fill", col.regions = bluered(256), xlab = "",
                  ylab = "", colorkey = colorkey)
    return(list(p = p, out.cn = as.data.frame(dat$out.cn),
                ploidy = ploidy, purity = purity, FGA = dat$FGA))
  }
  if (!is.null(seg)) {
    if (!is.null(outcome))
      names(outcome) <- samples
    if (!is.null(ordered))
      names(ordered) <- samples
    dat <- facets_dat(seg, filenames, path, samples, min_purity,
                      epsilon, adaptive)
    reducedM <- dat$out.cn
    samples <- samples[match(rownames(reducedM), samples)]
    if (!is.null(outcome))
      outcome <- outcome[match(rownames(reducedM), names(outcome))]
    if (!is.null(ordered) && !is.null(outcome))
      ordered <- order(outcome)
    rownames(reducedM) <- abbreviate(rownames(reducedM),
                                     minlength = 10)
    imagedata = reducedM
    imagedata[imagedata > 1.5] = 1.5
    imagedata[imagedata < -1.5] = -1.5
    if (is.null(ordered)) {
      cl = stats::hclust(stats::dist(imagedata), method = "ward")
      imagedata.ordered = imagedata[cl$order, ]
      imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
    if (!is.null(ordered)) {
      imagedata.ordered = imagedata[ordered, ]
      imagedata.ordered = as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
    chr = strsplit(colnames(imagedata), "\\.")
    chr = unlist(lapply(1:length(chr), function(x) chr[[x]][1]))
    chr = gsub("chr", "", chr)
    chr = as.numeric(chr)
    len = length(chr)
    chrom.ends <- rep(NA, length(table(chr)))
    d = 1
    for (r in unique(chr)) {
      chrom.ends[d] <- max(which(chr == r))
      d = d + 1
    }
    chrom.starts <- c(1, chrom.ends[-length(table(chr))] +
                        1)
    chrom.mids <- (chrom.starts + chrom.ends)/2
    bw = colorpanel(2, low = "white", high = "cadetblue4")
    colorkey = list(space = "right", height = 0.3, tick.number = 5)
    n <- nrow(reducedM)
    if (!is.null(outcome))
      x.lab <- outcome
    if (is.null(outcome))
      x.lab <- rep(" ", n)
    if (is.null(ordered))
      x.lab <- as.character(x.lab[cl$order])
    if (!is.null(ordered))
      x.lab <- as.character(x.lab[ordered])
    if (any(grepl("tcn", colnames(reducedM))) && any(grepl("ploidy",
                                                           colnames(reducedM))))
      scales = list(x = list(at = 1:n, labels = ploidy[cl$order],
                             rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                    z = list(at = n:1, labels = purity[cl$order],
                             rot = 90))
    else scales = list(x = list(at = 1:n, labels = x.lab,
                                rot = 90), y = list(at = len - chrom.mids, labels = names(table(chr))),
                       z = list(at = n:1, labels = rep(1, n), rot = 90))
    my.panel.levelplot.2 <- function(...) {
      panel.levelplot(...)
      panel.abline(h = len - chrom.starts[-1], col = "gray",
                   lwd = 1)
      panel.scales = list(x = list(at = 1:n), y = list(at = len -
                                                         chrom.mids), z = list())
    }
    my.panel = my.panel.levelplot.2
    p = levelplot(imagedata.ordered, panel = my.panel, scales = scales,
                  aspect = "fill", col.regions = bluered(256), xlab = "",
                  ylab = "", colorkey = colorkey)
    return(list(p = p, out.cn = as.data.frame(dat$out.cn)))
  }
}
