#' facets.heatmap
#'
#' Creates a heatmap of copy number profiles from segment files.
#'
#' @param seg a segmentation file containing the segmentation information of multiple patients
#' @param filenames the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param path the relative path to the files folder from your current directory
#' @param patients the names of the patients of the respective filenames. Default simply 1 to number of files.
#' @param min.purity the minimum purity of the sample required to be kept in the final dataset. Default is 0.3.
#' @param epsilon level of unions when aggregating segments between
#' @param ordered order in which patients should be printed. Default NUll leads to hierarchical clustering.
#' @param outcome for seg file only, if outcome associated with study it will be printed along the x axis for each patient
#' @return p a heatmap corresponding to the segnment files inputted
#' @export
#'
#' @examples library(gnomeR)
#'
#' @import
#' gplots
#' lattice


facets.heatmap <- function(seg = NULL,filenames = NULL, path, patients=NULL, min.purity = 0.3,
                           epsilon = 0.005,ordered = NULL,outcome = NULL){

  if(is.null(seg) && is.null(filenames))
    stop("You must provide either a complete segmentation file
         or a list of files to be loaded with their corresponding path")

  if(!is.null(seg) && !is.null(filenames))
    stop("Please provide either a complete segmentation file or a
         list of segmentation files to be loaded")

  ######################################

  if(!is.null(filenames)){
    dat <- facets.dat(seg = NULL,filenames, path, patients, min.purity, epsilon)
    reducedM <- dat$out.cn
    ploidy <- dat$ploidy
    purity <- dat$purity
    rownames(reducedM) <- abbreviate(rownames(reducedM),minlength = 10)
    imagedata=reducedM
    imagedata[imagedata>1.5]=1.5
    imagedata[imagedata< -1.5]= -1.5

    cl=hclust(dist(imagedata), method="ward")

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


  ###############################################################################

  if(!is.null(seg)){
    dat <- facets.dat(seg,patients = patients)
    reducedM <- dat$out.cn
    patients <- patients[match(rownames(reducedM),patients)]
    if(!is.null(outcome)) outcome <- outcome[match(rownames(reducedM),names(outcome))]
    if(!is.null(ordered) && !is.null(outcome)) ordered <- order(outcome)
    # if(!is.null(ordered) && !is.null(outcome)) {
    #   if(length(-which(is.na(match(patients,rownames(reducedM))))) > 0){
    #     outcome <- outcome[-which(is.na(match(patients,rownames(reducedM))))]
    #     ordered <- order(outcome)
    #     outcome <- outcome[ordered]
    #   }
    # }
    rownames(reducedM) <- abbreviate(rownames(reducedM),minlength = 10)
    imagedata=reducedM
    imagedata[imagedata>1.5]=1.5
    imagedata[imagedata< -1.5]= -1.5

    if(is.null(ordered)){
      cl=hclust(dist(imagedata), method="ward")

      imagedata.ordered=imagedata[cl$order,]
      imagedata.ordered=as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
    if(!is.null(ordered)){
      imagedata.ordered=imagedata[ordered,]
      imagedata.ordered=as.matrix(rev(as.data.frame(imagedata.ordered)))
    }
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

    n <- nrow(reducedM)

    if(!is.null(outcome)) x.lab <- outcome[ordered]
    else x.lab <- rep(1,n)
    # if(!is.null(ordered)) x.lab <- x.lab[ordered]
    if(is.null(ordered)) x.lab <- x.lab[cl$order]

    scales = list(x = list(at=1:n,labels=x.lab,rot=90),
                  y = list(at = len - chrom.mids, labels = names(table(chr))),
                  z = list(at=n:1,labels=rep(1,n),rot=90)) #[cl$order]


    my.panel.levelplot.2 <- function(...) {
      panel.levelplot(...)
      panel.abline(h = len - chrom.starts[-1], col = "gray", lwd = 1)
      panel.scales = list(x = list(at=1:n), y = list(at = len - chrom.mids), z = list())
    }
    my.panel = my.panel.levelplot.2

    p=levelplot(imagedata.ordered, panel = my.panel, scales=scales,aspect="fill",
                col.regions = bluered(256), xlab = "", ylab = "",colorkey=colorkey)

    return(list("p"=p,"out.cn"=dat$out.cn))
  }

}

