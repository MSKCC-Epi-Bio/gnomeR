#' facets.heatmap
#'
#' Creates a heatmap of copy number profiles from segment files.
#'
#' @param filenames the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param path the relative path to the files folder from your current directory
#' @param patients the names of the patients of the respective filenames. Default simply 1 to number of files.
#' @param min.purity the minimum purity of the sample required to be kept in the final dataset. Default is 0.3.
#' @param epsilon level of unions when aggregating segments between
#' @param ordered order in which patients should be printed. Default NUll leads to hierarchical clustering.
#' @return p a heatmap corresponding to the segnment files inputted
#' @export
#'
#' @examples library(gnomeR)
#'
#' @import
#' gplots
#' lattice


facets.heatmap <- function(filenames, path, patients=NULL, min.purity = 0.3,
                           epsilon = 0.005,ordered = NULL){

  dat <- facets.dat(filenames, path, patients, min.purity, epsilon)
  reducedM <- dat$out.cn
  ploidy <- dat$ploidy
  purity <- dat$purity
  rownames(reducedM) <- abbreviate(rownames(reducedM),minlength = 10)
  imagedata=reducedM
  imagedata[imagedata>1.5]=1.5
  imagedata[imagedata< -1.5]= -1.5

  cl=hclust(dist(imagedata), method="ward")
  plot(cl,hang= -1)

  imagedata.ordered=imagedata[cl$order,]
  imagedata.ordered=as.matrix(rev(as.data.frame(imagedata.ordered)))

  chr=strsplit(colnames(imagedata),"\\.")
  chr=unlist(lapply(1:length(chr),function(x)chr[[x]][1]))
  chr=gsub("chr","",chr)
  chr=as.numeric(chr)

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

  bw=colorpanel(2,low="white",high="cadetblue4")

  colorkey = list(space = "right", height = 0.3, tick.number = 5)

  n <- length(ploidy)
  #scales = list(x = list(at=1:n,labels=ploidy[cl$order],rot=90), y = list(at = len - chrom.mids, labels = names(table(chr))), z = list(at=n:1,labels=purity[cl$order],rot=90))
  scales = list(x = list(at=1:n,labels=ploidy[cl$order],rot=90),
                y = list(at = len - chrom.mids, labels = names(table(chr))),
                z = list(at=n:1,labels=purity[cl$order],rot=90))


  my.panel.levelplot.2 <- function(...) {
    panel.levelplot(...)
    #panel.abline(v = (cluster.start[-1] - 0.5), col = "black", lwd = 1, lty = 1)
    panel.abline(h = len - chrom.starts[-1], col = "gray", lwd = 1)
    panel.scales = list(x = list(at=1:n), y = list(at = len - chrom.mids), z = list())
  }
  my.panel = my.panel.levelplot.2

  p=levelplot(imagedata.ordered, panel = my.panel, scales=scales,
              col.regions = bluered(256), xlab = "", ylab = "",colorkey=colorkey)
  return(list("p"=p,"out.cn"=dat$out.cn,"ploidy"=ploidy,"purity"=purity,"FGA"=dat$FGA))
}

