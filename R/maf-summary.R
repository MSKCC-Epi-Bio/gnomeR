#' maf.summary
#' Creates a set of plot summarising a maf file.
#' @param maf the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param ... any argument belonging to the binmat method
#' @return p.class Barplot of counts of each variant classification
#' @return p.type Barplot of counts of each variant type
#' @return p.SNV Histogram of counts of each SNV class
#' @return p.patient.variant Histogram of counts variants per patient
#' @return p.variant.bp Boxplot of the distribution of variant classification per patient
#' @return p.variant.dist Boxplot of the relative proportion of each variant class in individual patients
#' @return p.variant.dist.bar Stacked barplot of the variants distribution in all patients
#' @return p.SNV.dist Boxplot of the relative proportion of each SNV class in individual patients
#' @return p.corr Correlation heatmap of the top 10 genes
#' @return p.comut Heatmap of the comutation of the top 10 genes
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' all.plots <- maf.summary(maf=mut %>% filter(Tumor_Sample_Barcode %in% patients) )
#' all.plots <- maf.summary(maf=mut %>%
#' filter(Tumor_Sample_Barcode %in% patients),spe.plat = TRUE)
#' @import
#' dplyr
#' dtplyr
#' stringr
#' ggplot2
#' reshape2
#' RColorBrewer
#' GGally


maf_viz <- function(maf,...){

  all_plots <- list(
    varclass = ggvarclass,
       vartype = ggvartype,
       nvclass = ggsnvclass,
      samplevar = ggsamplevar,
      topgenes = ggtopgenes) %>%
  invoke_map(, maf)

  return(all_plots)
}

ggvarclass <- function(maf) {

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


ggvartype <- function(maf) {

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


ggsnvclass <- function(maf) {

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


ggsamplevar <- function(maf) {

  maf2 <- maf %>%
   group_by(Tumor_Sample_Barcode) %>%
    mutate(n_alts = n()) %>%
    ungroup() %>%
    mutate(Tumor_Sample_Barcode = .data$Tumor_Sample_Barcode %>%
        forcats::fct_infreq())

  # distribution of variant per sample
  p.patient.variant <- maf2 %>%
    ggplot(aes(x = .data$Tumor_Sample_Barcode, fill = .data$Variant_Classification)) +
    geom_bar(position = "stack") +
    ggtitle("Variants per sample") +
    ylab("Variant Count") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  p.patient.variant
}


ggtopgenes <- function(maf, n_genes = 10) {

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
        forcats::fct_infreq() %>% fct_rev())

  p.genes <-  maf2 %>%
    ggplot(aes(x = .data$Hugo_Symbol,
               fill = .data$Variant_Classification)) +
    geom_bar(position = "stack") +
    coord_flip() +
    ggtitle("Top genes variants classification") + xlab("Gene Name")

  p.genes
}

ggvarprop <- function(mut) {

  nb.cols = 20

  varprop.maf <- calc_variant_prop(mut)

  p.variant.dist <- varprop.maf %>%
    ggplot(aes(x = .data$Variant_Classification,y = .data$varProp)) +
    geom_boxplot(aes(fill = .data$Variant_Classification)) +
    scale_fill_manual(values =
                        grDevices::colorRampPalette(
                          RColorBrewer::brewer.pal(8, "Accent"))(nb.cols)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none") + ylab("% Variant")

  p.variant.dist

}

ggvarprop2 <- function(mut) {

  nb.cols = 20

  varprop.maf <- calc_variant_prop(mut)

  p.variant.dist.bar <- varprop.maf %>%
    ggplot(aes(x = .data$Tumor_Sample_Barcode, y=.data$varProp)) +
    geom_bar(aes(fill = .data$Variant_Classification),stat = "identity") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_fill_manual(values = grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(8, "Accent"))(nb.cols) ) +
    ylab("% Variants")

  p.variant.dist.bar

}

ggsnvclass <- function(mut) {

  nb.cols = 20

  maf <- maf %>%
    filter(
      .data$Variant_Type == "SNP",
      .data$HGVSc != ""
    ) %>%
    mutate(SNV_Class = substrRight(.data$HGVSc, 3)) %>%
    mutate(SNV_Class = .data$SNV_Class %>%
        forcats::fct_infreq() %>%
        forcats::fct_rev())


    p.SNV.dist <- maf %>%
      group_by(.data$Tumor_Sample_Barcode) %>%
      mutate(totalMut = n()) %>%
      ungroup() %>%
      group_by(.data$Tumor_Sample_Barcode, .data$SNV_Class) %>%
      summarise(N=n(),
                varProp = .data$N/unique(.data$totalMut)) %>%
      ggplot(aes(x = .data$SNV_Class,y = .data$varProp)) + geom_boxplot(aes(fill = .data$SNV_Class))+
      scale_fill_manual(values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Accent"))(nb.cols)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position="none") + ylab("% SNV")

    p.SNV.dist
}

gggenecor <- function(mut, ...) {

  bin.maf <- binmat(maf = maf,...)

  keep <- names(sort(apply(bin.maf,2,function(x){sum(x)}),decreasing = T))[1:10]
  bin.maf <- bin.maf[,keep]

  p.corr <- GGally::ggcorr(bin.maf,limits = NULL)

  p.corr

}

ggcomut <- function(mut, ...) {

  bin.maf <- binmat(maf = maf,...)
  co.mut <- apply(bin.maf,2,function(x){
    apply(bin.maf,2,function(y){
      sum(y == 1 & x == 1,na.rm = T)/length(x)
    })
  })

  p.comut <- GGally::ggcorr(co.mut,limits = NULL)

  p.comut
}





