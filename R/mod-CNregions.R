#' CNregions.mod
#' Modified CNregions function from the facets package to handle small sample sizes.
#' @param seg a segmentation file containing the segmentation information of multiple patients
#' @param epsilon the maximum Euclidean distance between adjacent probes tolerated for denying a nonredundant region. epsilon=0 is equivalent to taking the union of all unique break points across the n samples.
#' @param adaptive Vector of length-m lasso penalty terms.
#' @param rmCNV If TRUE, remove germline CNV.
#' @param cnv A data frame containing germline CNV data.
#' @param frac.overlap A parameter needed to be explain.
#' @param rmSmallseg If TRUE, remove small segment.
#' @param nProbes The segment length threshold below which the segment will be removed if rmSmallseq = TRUE.
#' @return reducedM : A matrix ready for plotting.
#' @export
#' @examples library(gnomeR)
#' CNregions.mod(seg[1:50])
CNregions.mod <- function (seg, epsilon = 0.005, adaptive = FALSE, rmCNV = FALSE,
                           cnv = NULL, frac.overlap = 0.5, rmSmallseg = TRUE, nProbes = 15)
{
  colnames(seg) = c("sample", "chromosome", "start", "end",
                    "num.mark", "seg.mean")
  seg <- seg %>% mutate(chromosome = as.numeric(.data$chromosome))
  seg = subset(seg, seg$chromosome <= 22)
  seg = seg[order(seg[, 1], seg[, 2], seg[, 3]), ]
  if (rmSmallseg) {
    seg = subset(seg, seg$num.mark >= nProbes)
  }
  if (rmCNV) {
    cat("Removing CNV...", "\n")
    if (is.null(cnv))
      stop("To remove CNV, please include the cnv file. Otherwise set rmCNV=F")
    gr.cnv = GenomicRanges::GRanges(seqnames = as.character(cnv[, 1]), ranges = IRanges::IRanges(start = as.numeric(cnv[,
                                                                                                2]), end = as.numeric(cnv[, 3])))
    gr.seg = GenomicRanges::GRanges(seqnames = paste("chr", seg[, 2], sep = ""),
                     ranges = IRanges::IRanges(start = seg[, 3], end = seg[, 4]))
    overlap = GenomicRanges::findOverlaps(gr.cnv, gr.seg)
    queryHits = queryHits(overlap)
    subjectHits = subjectHits(overlap)
    start = pmax(seg[subjectHits, 3], cnv[queryHits, 2])
    end = pmin(seg[subjectHits, 4], cnv[queryHits, 3])
    seg.length = seg[subjectHits, 4] - seg[subjectHits, 3]
    seg.length[seg.length == 0] = 1
    p.overlap = (end - start)/seg.length
    cnv.idx = subjectHits[which(p.overlap >= frac.overlap)]
    cnv.idx = unique(cnv.idx)
    seg = seg[-cnv.idx, ]
  }
  if (dim(seg)[1] == 0)
    stop("seg file is empty.")
  samples = unique(seg$sample)
  chr = start.maploc = end.maploc = nur = out = NULL
  for(i in unique(as.numeric(seg[, 2]))) { #unique(as.numeric(seg[, 2]))
    subdata = subset(seg, seg$chromosome == i)
    u.breakpts = sort(unique(subdata[, 3]))
    l = length(u.breakpts)
    outi = NULL
    for (j in 1:length(samples)) {
      sj = subset(subdata, sample == samples[j])
      if (dim(sj)[1] == 0)
        stop("Found sample(s) with no segments. Check parameter setting. Relax threshold values.")
      outj = rep(NA, l)
      sj$start[1] = u.breakpts[1]
      idx = c(match(sj$start, u.breakpts), l)
      outj = c(rep(sj[, 6], times = diff(idx)), utils::tail(sj[,
                                                        6], 1))
      outi = cbind(outi, outj)
    }
    # print(paste0(i,",",dim(outi)))
    u.start = u.breakpts
    u.se = sort(union(unique(subdata$end), unique(subdata$start)))
    u.end = lapply(2:length(u.start), FUN = function(x) u.se[which(u.start[x] <=
                                                                     u.se)[1]])
    u.end = unlist(u.end)
    u.end = u.end[stats::complete.cases(u.end)]
    u.end = c(u.end, max(u.se))
    if(nrow(outi) > 1) adjrow.dist = apply(sqrt(diff(outi)^2), 1, mean)
    else adjrow.dist = mean(sqrt(diff(outi[1,])^2))
    adjrow.dist = c(0, adjrow.dist)
    if (adaptive) {
      epsilon = stats::quantile(adjrow.dist, prob = 0.97)
    }
    seg.break = which(adjrow.dist >= epsilon)
    if (length(seg.break) == 0)
      stop("Check parameter setting. Set lower epsilon value.")
    start = seg.break
    if(length(u.start) > 1){
      if (start[1] > 1) {
        start = c(1, start)
      }
      end = seg.break - 1
      if (end[1] == 0) {
        end = end[-1]
      }
      len = dim(outi)[1]
      if (tail(end, 1) < len) {
        end = c(end, len)
      }
    }
    else{
      start = 1
      end = 1
    }
    nuri = end - start + 1
    rowidx = rep(c(1:length(nuri)), nuri)
    nur = c(nur, nuri)
    get.medoid = function(x) {
      if (dim(x)[1] > 1) {
        pam(x, k = 1)$medoid
      }
      else {
        x
      }
    }
    outi.medoid = by(outi, rowidx, get.medoid)
    outi.medoid = matrix(unlist(outi.medoid), nrow = length(outi.medoid),
                         byrow = T)
    start.maploc = c(start.maploc, u.start[start])
    end.maploc = c(end.maploc, u.end[end])
    chr = c(chr, rep(i, dim(outi.medoid)[1]))
    out = rbind(out, outi.medoid)
  }
  colnames(out) = samples
  rownames(out) = paste("chr", chr, ".", start.maploc, "-",
                        end.maploc, sep = "")
  reducedM = t(out)
  return(reducedM)
}
