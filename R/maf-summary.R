#' maf.summary
#'
#' Creates a set of plot summarising a maf file.
#'
#' @param maf the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @return A series of plot summarising information about the maf file inputted.
#' @export
#'
#' @examples library(gnomeR)
#' all.plots <- maf.summary(maf=mut)
#' @import
#' dplyr
#' stringr
#' ggplot2
#' reshape2


maf.summary <- function(maf,mut.type = "SOMATIC"){

  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
  if(length(match("Variant_Classification",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
  if(length(match("Mutation_Status",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a mutation status column. (Mutation_Status)")

  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
  if (sum(grepl("KMT2D", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "KMT2D" ~ "MLL2",
        TRUE ~ Hugo_Symbol
      ))

    warning("KMT2D has been recoded to MLL2")
  }

  if (sum(grepl("KMT2C", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "KMT2C" ~ "MLL3",
        TRUE ~ Hugo_Symbol
      ))

    warning("KMT2C has been recoded to MLL3")
  }

  substrRight <- function(x, n){
    x <- as.character(x)
    substr(x, nchar(x)-n+1, nchar(x))
  }

  # maf <- as.data.frame(maf)
  if(mut.type == "ALL") Mut.filt = unique(maf$Mutation_Status) else Mut.filt = mut.type
  maf <- maf %>% filter(tolower(Mutation_Status) %in% tolower(Mut.filt))

  ## summarise variant wise ##

  # variants call summary plot #
  maf$Variant_Classification <- factor(as.character(maf$Variant_Classification),
                                       levels =  as.character(unlist(maf %>%
                                                                       group_by(Variant_Classification) %>%
                                                                       summarise(N = n()) %>%
                                                                       arrange(-N) %>% select(Variant_Classification))))

  p.class <- maf %>% ggplot(aes(x = Variant_Classification,color=Variant_Classification,fill = Variant_Classification)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("Variant Classification count") + xlab("Variant Classification")


  # variants type summary plot #
  maf$Variant_Type <- factor(as.character(maf$Variant_Type),
                             levels =  as.character(unlist(maf %>%
                                                             group_by(Variant_Type) %>%
                                                             summarise(N = n()) %>%
                                                             arrange(-N) %>% select(Variant_Type))))
  p.type <- maf %>% ggplot(aes(x = Variant_Type,color=Variant_Type,fill = Variant_Type)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("Variant Type count") + xlab("Variant Type")

  # SNV #
  p.SNV <- maf %>%
    filter(Variant_Type == "SNP",
           HGVSc != "") %>%
    mutate(SNV_Class = substrRight(HGVSc,3)) %>%
    filter(!grepl("N",SNV_Class)) %>%
    ggplot(aes(x = SNV_Class,color=SNV_Class,fill = SNV_Class)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("SNV Class count") + xlab("SNV Class")



  ## summarise patient wise ##
  # distribution of variant per sample #
  maf$Tumor_Sample_Barcode <- factor(as.character(maf$Tumor_Sample_Barcode),levels = as.character(unlist(maf %>%
                                                                                                           group_by(Tumor_Sample_Barcode) %>%
                                                                                                           summarise(N = n()) %>%
                                                                                                           arrange(-N) %>% select(Tumor_Sample_Barcode))))

  p.patient.variant <- maf %>%
    group_by(Tumor_Sample_Barcode,Variant_Classification) %>%
    ggplot(aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
    ggtitle("Variants per sample") + ylab("Variant Count") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  # same but in boxplot format # --> simpler to read
  p.variant.bp <- maf %>%
    group_by(Tumor_Sample_Barcode,Variant_Classification) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Variant_Classification, y = N,color = Variant_Classification)) +
    ggtitle("Variant classification distribution") + xlab("Variant Classification") +
    geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


  # most mutated genes and class distrib #
  top.genes <- as.character(unlist(maf %>%
                                     group_by(Hugo_Symbol) %>%
                                     summarise(N = n()) %>%
                                     arrange(-N) %>%
                                     select(Hugo_Symbol)))
  maf$Hugo_Symbol <- factor(as.character(maf$Hugo_Symbol),levels = top.genes)

  p.genes <-  maf %>%
    filter(Hugo_Symbol %in% top.genes[1:10]) %>%
    ggplot(aes(x = Hugo_Symbol)) + geom_bar(aes(fill = Variant_Classification)) +
    coord_flip() +
    ggtitle("Top genes variants classification") + xlab("Gene Name")


  return(list("p.class"=p.class,"p.type"=p.type,"p.SNV" = p.SNV,
    "p.patient.variant"=p.patient.variant,"p.variant.bp" = p.variant.bp,"p.genes" = p.genes))
}
