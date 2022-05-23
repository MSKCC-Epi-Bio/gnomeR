#' Creates a set of plot summarising a mutation file.
#' @param mutation Raw mutation dataframe containing alteration data
#' @param ... any argument belonging to the gene_binary method
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
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' all.plots <- mutation_viz(mutation=mut %>% filter(Tumor_Sample_Barcode %in% samples))
#' all.plots <- mutation_viz(mutation=mut %>%
#' filter(Tumor_Sample_Barcode %in% samples),specify_panel = TRUE)
#' @import
#' dplyr
#' dtplyr
#' stringr
#' ggplot2
#' GGally
#' ComplexHeatmap

mutation_viz <- function(mutation, ...) {


    all_plots <- list(
      varclass = ggvarclass,
      vartype = ggvartype,
      snvclass = ggsnvclass,
      samplevar = ggsamplevar,
      topgenes = ggtopgenes,
      genecor = gggenecor) %>%
      purrr::invoke_map(, mutation)

  return(all_plots)
}



#' Add a percentage to counts
#'
#' @param x A barplot of the mutationviz plot that follow barplot visualization
#' @param ... other arguments as passed to adjust the percentage label size
#' @return mutationviz Barplot The same barplot is now returned with percentages
#' @export
#'
#' @examples
#' ggvarclass(mut) + add.perc()
#'

add.perc<-function(x,...){geom_text(
  aes(label=paste0(round(stat(.data$prop)*100,1),"%"), group=1),
  stat="count",
  hjust=0, nudge_y = -0.25,...)}


#' Barplot of Variant Classification Counts
#'
#' @param mutation Raw mutation dataframe containing alteration data
#'
#' @return Barplot of counts of each variant classification
#' @export
#'
#' @examples
#' ggvarclass(mut)
#'
ggvarclass <- function(mutation) {


  # relevel Variant Classification by frequency
  mutation <- mutation %>%
    mutate(Variant_Classification =
             stringr::str_replace_all(.data$Variant_Classification, "_", " ")) %>%
    mutate(Variant_Classification = .data$Variant_Classification %>%
             forcats::fct_infreq() %>%
             forcats::fct_rev())

  p.class <- mutation %>%
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
#' @param mutation Raw mutation dataframe containing alteration data
#'
#' @return Barplot of counts of each variant type
#' @export
#'
#' @examples
#' ggvartype(mut)
#'
ggvartype <- function(mutation) {

  # relevel Variant Type by frequency
  mutation <- mutation %>%
    mutate(Variant_Type = .data$Variant_Type %>%
             forcats::fct_infreq() %>%
             forcats::fct_rev())

  p.type <- mutation %>%
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
#' @param mutation Raw mutation dataframe containing alteration data
#'
#' @return Histogram of counts of each SNV class
#' @export
#'
#' @examples
#' ggsnvclass(mut)
#'
ggsnvclass <- function(mutation) {


  # filter only SNPs
  mutation <- mutation %>%
    filter(
      .data$Variant_Type == "SNP",
      .data$HGVSc != ""
    ) %>%
    mutate(SNV_Class = substrRight(.data$HGVSc, 3)) %>%
    mutate(SNV_Class = .data$SNV_Class %>%
             forcats::fct_infreq() %>%
             forcats::fct_rev())

  p.SNV <- mutation %>%
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
#' @param mutation Raw mutation dataframe containing alteration data
#'
#' @return Histogram of counts of variants per tumor sample
#' @export
#'
#' @examples
#' ggsamplevar(mut)
#'
ggsamplevar <- function(mutation) {

  mutation2 <- mutation %>%
    group_by(.data$Tumor_Sample_Barcode) %>%
    mutate(n_alts = n()) %>%
    ungroup() %>%
    mutate(Tumor_Sample_Barcode = .data$Tumor_Sample_Barcode %>%
             forcats::fct_infreq())

  # distribution of variant per sample
  p.patient.variant <- mutation2 %>%
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
#' @param mutation Raw mutation dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @return Barplot of counts of top variant genes
#' @export
#'
#' @examples
#' ggtopgenes(mut)
#'
ggtopgenes <- function(mutation, n_genes = 10) {


  top_genes <- mutation %>%
    group_by(.data$Hugo_Symbol) %>%
    summarise(N = n()) %>%
    arrange(-.data$N) %>%
    select(.data$Hugo_Symbol) %>%
    pull(.data$Hugo_Symbol)

  top_genes <- top_genes[1:min(length(top_genes),n_genes)] %>%
    as.character()

  mutation2 <- mutation %>%
    filter(.data$Hugo_Symbol %in% top_genes) %>%
    ungroup() %>%
    mutate(Hugo_Symbol = .data$Hugo_Symbol %>%
             forcats::fct_drop() %>%
             forcats::fct_infreq() %>%
             forcats::fct_rev())

  p.genes <-  mutation2 %>%
    ggplot(aes(x = .data$Hugo_Symbol,
               fill = .data$Variant_Classification)) +
    geom_bar(position = "stack") +
    coord_flip() +
    ggtitle("Top genes variants classification") + xlab("Gene Name")

  p.genes
}

#' Correlation Heatmap of the Top Altered Genes
#'
#' @param mutation Raw mutation dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @param ... Further create_gene_binary() arguments
#' @return Correlation heatmap of the top altered genes
#' @export
#'
#' @examples
#' gggenecor(mut)
#'
gggenecor <- function(mutation, n_genes = 10, ...) {

  bin.mutation <- create_gene_binary(mutation = mutation,...)

  keep <- names(sort(apply(bin.mutation,2,
                           function(x){sum(x)}),
                     decreasing = T))
  keep <- keep[1:min(length(keep),n_genes)]
  bin.mutation <- bin.mutation[,keep]

  p.corr <- GGally::ggcorr(dat = bin.mutation, cor_matrix = stats::cor(bin.mutation),limits = NULL)

  p.corr

}





#'
#' Comutation Heatmap of the Top Altered Genes
#'
#' @param mutation Raw mutation dataframe containing alteration data
#' @param n_genes Number of top genes to display in plot
#' @param ... Further create_gene_binary() arguments
#' @return Comutation heatmap of the top genes
#' @export
#'
#' @examples
#' ggcomut(gnomeR::mut)
#'
ggcomut <- function(mutation, n_genes = 10, ...) {

  bin.mutation <- create_gene_binary(mutation = mutation,...)
  keep <- names(sort(apply(bin.mutation,2,
                           function(x){sum(x)}),
                     decreasing = T))
  keep <- keep[1:min(length(keep),n_genes)]
  bin.mutation <- bin.mutation[,keep]

  co.mut <- apply(bin.mutation,2,function(x){
    apply(bin.mutation,2,function(y){
      sum(y == 1 & x == 1,na.rm = T)/length(x)
    })
  })

  p.comut <- GGally::ggcorr(dat = bin.mutation, cor_matrix = co.mut, limits = NULL)

  p.comut
}

#' Heatmap of all events after gene_binary - using binary distance
#'
#' @param hmat dataframe obtained after create_gene_binary()
#' @param ... Further arguments as passed to ComplexHeatmap::Heatmap
#' @return heatmap of gene_binary events
#' @export
#'
#' @examples
#' set.seed(123)
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' gen_dat <- create_gene_binary(samples=samples, mutation=mut, cna=cna, fusion=fusion, set.plat=TRUE)
#' ggheatmap(gen_dat, show_row_names=FALSE, show_column_names=FALSE)
#'

ggheatmap <- function(hmat, ...){

  #check if the matrix is not binary
  if(sum(hmat==0, na.rm=T) +
     sum(hmat==1, na.rm=T) +
     sum(is.na(hmat)) != (nrow(hmat) * ncol(hmat))) {
    stop("ggheatmap can only be plotted when gene_binary is binary, set cna.binary=TRUE")
  }

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


  hmat.colors = structure(c("white","black", "coral","cadetblue", "forestgreen"),
                          names=c("0","1","2","3","4"))

  hmap.legend = list(
    title = "Events",
    at = c(0,1,2,3,4),
    labels = c("WT", "Mut", "Amplification","Deletion", "Fusion")
  )

  ComplexHeatmap::Heatmap(t(hmat),col=hmat.colors , na_col="grey",
                          clustering_distance_rows="binary",clustering_distance_columns="binary", heatmap_legend_param = hmap.legend, ...)

}



