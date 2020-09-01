
#' Creates a set of plot summarising a maf file.
#' @param maf Raw maf dataframe containing alteration data
#' @param ... any argument belonging to the binmat method
#' @return Returns a list of the following plots:
#' @return varclass Barplot of counts of each variant classification
#' @return vartype Barplot of counts of each variant type
#' @return snvclass Histogram of counts of each SNV class
#' @return samplevar Histogram of counts variants per patient
#' @return topgenes Barplot of counts of top variant genes
#' @return genecor Correlation heatmap of the top 10 genes
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' all.plots <- maf_viz(maf=mut %>% filter(Tumor_Sample_Barcode %in% patients))
#' all.plots <- maf_viz(maf=mut %>%
#' filter(Tumor_Sample_Barcode %in% patients),spe.plat = TRUE)
#' @import
#' dplyr
#' dtplyr
#' stringr
#' ggplot2
#' GGally
#' ComplexHeatmap

maf_viz <- function(maf, ...) {
  all_plots <- list(
    varclass = ggvarclass,
    vartype = ggvartype,
    snvclass = ggsnvclass,
    samplevar = ggsamplevar,
    topgenes = ggtopgenes,
    genecor = gggenecor) %>%
    purrr::invoke_map(, maf)

  return(all_plots)
}


add.perc<-function(x,...){geom_text(
       aes(label=paste0(round(stat(prop)*100,2),"%"), group=1),
       stat="count",
       hjust=0, nudge_y = -0.25,...)}


#' Barplot of Variant Classification Counts
#'
#' @param maf Raw maf dataframe containing alteration data
#'
#' @return Barplot of counts of each variant classification
#' @export
#'
#' @examples
#' ggvarclass(mut)
#'
ggvarclass <- function(maf) {

  check_maf_column(maf = maf, "Variant_Classification")

  # relevel Variant Classification by frequency
  maf <- maf %>%
    mutate(Variant_Classification =
        stringr::str_replace_all(.data$Variant_Classification, "_", " ")) %>%
    mutate(Variant_Classification = .data$Variant_Classification %>%
        forcats::fct_infreq() %>%
        forcats::fct_rev())

  p.class <- maf %>%
    ggplot(aes(x = .data$Variant_Classification,
      fill = .data$Variant_Classification)) +
    geom_bar() +
    coord_flip() +
    theme(legend.position="none") +
    ggtitle("Variant Classification Count") +
    xlab("Variant Classification")

  p.class
}

#' Barplot of Variant Type Counts
#'
#' @param maf Raw maf dataframe containing alteration data
#'
#' @return Barplot of counts of each variant type
#' @export
#'
#' @examples
#' ggvartype(mut)
#'
ggvartype <- function(maf) {

  check_maf_column(maf = maf, "Variant_Type")

  # relevel Variant Type by frequency
  maf <- maf %>%
    mutate(Variant_Type = .data$Variant_Type %>%
        forcats::fct_infreq() %>%
        forcats::fct_rev())

  p.type <- maf %>%
    ggplot(aes(x = .data$Variant_Type,
      color=.data$Variant_Type,
      fill = .data$Variant_Type)) +
    geom_bar() +
    coord_flip() +
    theme_minimal() +
    theme(legend.position="none") +
    ggtitle("Variant Type Count") +
    xlab("Variant Type")

  p.type

}

#' Histogram of SNV class Counts
#'
#' @param maf Raw maf dataframe containing alteration data
#'
#' @return Histogram of counts of each SNV class
#' @export
#'
#' @examples
#' ggsnvclass(mut)
#'
ggsnvclass <- function(maf) {

  # NOTE: should make function accept multiple args
  check_maf_column(maf = maf, "HGVSc")
  check_maf_column(maf = maf, "Variant_Type")

  # filter only SNPs
  maf <- maf %>%
    filter(
      .data$Variant_Type == "SNP",
      .data$HGVSc != ""
    ) %>%
    mutate(SNV_Class = substrRight(.data$HGVSc, 3)) %>%
    mutate(SNV_Class = .data$SNV_Class %>%
        forcats::fct_infreq() %>%
        forcats::fct_rev())

  p.SNV <- maf %>%
    filter(!grepl("N", .data$SNV_Class)) %>%
    ggplot(aes(x = .data$SNV_Class, color = .data$SNV_Class,
      fill = .data$SNV_Class)) +
    geom_bar() +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle("SNV Class Count") +
    xlab("SNV Class")

    p.SNV
}

#' Histogram of Variants Per Sample Colored By Variant Classification
#'
#' @param maf Raw maf dataframe containing alteration data
#'
#' @return Histogram of counts of variants per tumor sample
#' @export
#'
#' @examples
#' ggsamplevar(mut)
#'
ggsamplevar <- function(maf) {

  check_maf_column(maf = maf, "Variant_Classification")
  check_maf_column(maf = maf, "Tumor_Sample_Barcode")

  maf2 <- maf %>%
   group_by(Tumor_Sample_Barcode) %>%
    mutate(n_alts = n()) %>%
    ungroup() %>%
    mutate(Tumor_Sample_Barcode = .data$Tumor_Sample_Barcode %>%
        forcats::fct_infreq())

  # distribution of variant per sample
  p.patient.variant <- maf2 %>%
    ggplot(aes(x = .data$Tumor_Sample_Barcode,
               fill = .data$Variant_Classification)) +
    geom_bar(position = "stack") +
    ggtitle("Variants per sample") +
    ylab("Variant Count") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  p.patient.variant
}


#' Barplot of Most Frequently Altered Genes
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @return Barplot of counts of top variant genes
#' @export
#'
#' @examples
#' ggtopgenes(mut)
#'
ggtopgenes <- function(maf, n_genes = 10) {

  check_maf_column(maf = maf, "Hugo_Symbol")

  top_genes <- maf %>%
    group_by(.data$Hugo_Symbol) %>%
    summarise(N = n()) %>%
    arrange(-.data$N) %>%
    select(.data$Hugo_Symbol) %>%
    pull(Hugo_Symbol)

  top_genes <- top_genes[1:n_genes] %>%
    as.character()

    maf2 <- maf %>%
      filter(.data$Hugo_Symbol %in% top_genes) %>%
      ungroup() %>%
      mutate(Hugo_Symbol = .data$Hugo_Symbol %>%
        forcats::fct_drop() %>%
        forcats::fct_infreq() %>%
          forcats::fct_rev())

  p.genes <-  maf2 %>%
    ggplot(aes(x = .data$Hugo_Symbol,
               fill = .data$Variant_Classification)) +
    geom_bar(position = "stack") +
    coord_flip() +
    ggtitle("Top genes variants classification") + xlab("Gene Name")

  p.genes
}

#' Correlation Heatmap of the Top Altered Genes
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @param ... Further binmat() arguments
#' @return Correlation heatmap of the top altered genes
#' @export
#'
#' @examples
#' gggenecor(mut)
#'
gggenecor <- function(maf, n_genes = 10, ...) {

  bin.maf <- binmat(maf = maf,...)

  keep <- names(sort(apply(bin.maf,2,
                           function(x){sum(x)}),
                     decreasing = T))[1:n_genes]
  bin.maf <- bin.maf[,keep]

  p.corr <- GGally::ggcorr(dat = bin.maf, cor_matrix = cor(bin.maf),limits = NULL)

  p.corr

}

#' Co-mutation Heatmap of the Top Altered Genes
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @param ... Further binmat() arguments
#' @return Comutation heatmap of the top genes
#' @export
#'
#' @examples
#' ggcomut(mut)
#'
ggcomut <- function(maf, n_genes = 10, ...) {

  bin.maf <- binmat(maf = maf,...)
  keep <- names(sort(apply(bin.maf,2,
                           function(x){sum(x)}),
                     decreasing = T))[1:n_genes]
  bin.maf <- bin.maf[,keep]

  co.mut <- apply(bin.maf,2,function(x){
    apply(bin.maf,2,function(y){
      sum(y == 1 & x == 1,na.rm = T)/length(x)
    })
  })

  p.comut <- GGally::ggcorr(dat = bin.maf, cor_matrix = co.mut, limits = NULL)

  p.comut
}

#' Heatmap of all events after binmat
#'
#' @param hmat dataframe obtained after binamt
#' @param ... Further arguments as passed to ComplexHeatmap::Heatmap
#' @return heatmap of binmat events
#' @export
#'
#' @examples
#' set.seed(123)
#' patients<-as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 100, replace=FALSE)]
#' ggheatmap(binmat(patients=patients, maf=mut, cna=cna, fusion=fusion, set.plat=TRUE), show_row_names=FALSE, show_column_names=FALSE)
#'

ggheatmap<-function(hmat, ...){

  #check if the matrix is not binary
  if(sum(hmat==0, na.rm=T) + sum(hmat==1, na.rm=T) + sum(is.na(hmat)) != (nrow(hmat) * ncol(hmat))) {stop("ggheatmap can only be plotted when binmat is binary, set cna.binary=TRUE")}

  idx.amp = grep(".Amp", colnames(hmat))
  if(length(idx.amp)>0){
    tt<- hmat[,idx.amp]; tt[ tt==1 ] <- 2
    hmat[,idx.amp]<-tt
  }

  idx.del = grep(".Del", colnames(hmat))
  if(length(idx.del)>0){
    tt<- hmat[,idx.del]; tt[tt==1] <- 3
    hmat[,idx.del]<-tt
  }


  idx.fus = c(grep(".fus", colnames(hmat)))
  if(length(idx.fus)>0){
    tt<- hmat[,idx.fus]; tt[tt==1] <- 4
    hmat[,idx.fus]<-tt
  }


  hmat.colors = structure(c("white","black", "coral","cadetblue", "forestgreen"), names=c("0","1","2","3","4"))

  hmap.legend = list(
    title = "Events",
    at = c(0,1,2,3,4),
    labels = c("WT", "Mut", "Amplification","Deletion", "Fusion")
  )

  ComplexHeatmap::Heatmap(t(hmat),col=hmat.colors , na_col="grey", clustering_distance_rows="binary",clustering_distance_columns="binary", heatmap_legend_param = hmap.legend, ...)

}



