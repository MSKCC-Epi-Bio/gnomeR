#' create.bin.matrix
#'
#' Enables creation of a binary matrix from a maf file with
#' a predifined list of patients (rows are patients and columns are genes)
#'
#' @param patients a character vector that let's the user specify the patients to be used to create the matrix.
#' Default is NULL is which case all patients in the MAF file will be used.
#' @param maf A MAF file.
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param SNP.only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include.silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @param fusion An optional MAF file for fusions. If inputed the outcome will be added to the matrix with columns ending in ".fus".
#' Default is NULL.
#' @param cna An optional CNA files. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
#' Default is NULL. Note that this file must have patients as columns and genes as rows. create.bin.matrix expects a matrix with
#' values between -2 and 2. Please do not use any other format. Other functions in the package are available to deal with more detailed
#' CNA data.
#' @param cna.relax for cna data only enables to count both gains and shallow deletions as amplifications and deletions respectively.
#' @param spe.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is F
#' @return mut : a binary matrix of mutation data
#' @return no.mu.patients : a character vector of patients having no mutations found in the MAF file.
#' @export
#' @examples library(gnomeR)
#' mut.only <- create.bin.matrix(maf = mut)
#' all.platforms <- create.bin.matrix(patients = unique(mut$Tumor_Sample_Barcode)[1:100],maf = mut,fusion = fusion,cna = cna)
#' @import dplyr
#' @import stringr

create.bin.matrix <- function(patients=NULL, maf, mut.type = "SOMATIC",SNP.only = F,include.silent = F,
                              fusion = NULL,cna = NULL,cna.relax = F, spe.plat = F){

  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
  if(length(match("Variant_Classification",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
  if(length(match("Mutation_Status",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a mutation status column. (Mutation_Status)")
  maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
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

  maf <- as.data.frame(maf)

  # filter/define patients #
  if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
  else patients <- as.character(unique(maf$Tumor_Sample_Barcode))

  # clean gen dat #
  if(SNP.only) SNP.filt = "SNP" else SNP.filt = unique(maf$Variant_Type)
  if(!include.silent) Variant.filt = "Silent" else Variant.filt = ""
  if(tolower(mut.type) == "all") Mut.filt = unique(maf$Mutation_Status) else Mut.filt = mut.type

  maf <- maf %>% filter(Variant_Classification != Variant.filt,
                        Variant_Type %in% SNP.filt,
                        tolower(Mutation_Status) %in% tolower(Mut.filt))


  #### out frame
  mut <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(maf$Hugo_Symbol))))
  colnames(mut) <- unique(maf$Hugo_Symbol)
  rownames(mut) <- patients

  for(i in patients){
    genes <- maf$Hugo_Symbol[maf$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0){mut[match(i,rownames(mut)),match(unique(as.character(genes)),colnames(mut))] <- 1}
  }

  missing.mut <- apply(mut,1,function(x){sum(x)==0})
  if(sum(missing.mut) > 0)
    warning("Some patients did not have any mutations found in the MAF file.")


  # add fusions if needed #
  if(!is.null(fusion)){
    fusion <- as.data.frame(fusion)
    # quick data checks #
    if(length(match("Tumor_Sample_Barcode",colnames(fusion))) == 0)
      stop("The fusion file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
    if(length(match("Hugo_Symbol",colnames(fusion))) == 0)
      stop("The fusion file inputted is missing a gene name column. (Hugo_Symbol)")


    fusion <- fusion %>%
      filter(Tumor_Sample_Barcode %in% patients)

    #### out frame
    fusion.out <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(fusion$Hugo_Symbol))))
    colnames(fusion.out) <- unique(fusion$Hugo_Symbol)
    rownames(fusion.out) <- patients

    for(i in patients){
      genes <- fusion$Hugo_Symbol[fusion$Tumor_Sample_Barcode %in% i]
      if(length(genes) != 0){fusion.out[match(i,rownames(fusion.out)),
                                        match(unique(as.character(genes)),colnames(fusion.out))] <- 1}
    }
    colnames(fusion.out) <- paste0(colnames(fusion.out),".fus")
    mut <- as.data.frame(cbind(mut,fusion.out))
    rownames(mut) <- patients
  }


  # add CNA if needed #
  if(!is.null(cna)){
    cna <- as.data.frame(cna)
    rownames(cna) <- cna[,1]
    cna <- cna[,-1]
    cna <- as.data.frame(t(cna))
    rownames(cna) <- gsub("\\.","-",rownames(cna))
    cna <- cna[rownames(cna) %in% patients,]

    temp <- do.call("cbind",apply(cna,2,function(x){
      if(cna.relax){
        yA <- ifelse(x>=0.9,1,0)
        yD <- ifelse(x<=-0.9,1,0)
      }
      if(!cna.relax){
        yA <- ifelse(x==2,1,0)
        yD <- ifelse(x==-2,1,0)
      }
      out <- as.data.frame(cbind(yA,yD))
      colnames(out) <- c("Amp","Del")
      return(out)
    }))

    cna <- temp[,apply(temp,2,function(x){sum(x,na.rm=T) > 0})]

    # add missing
    missing <- patients[which(is.na(match(patients,rownames(cna))))]
    add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
    rownames(add )  <- missing
    colnames(add )<- colnames(cna)
    cna <- as.data.frame(rbind(cna,add))
    cna <- cna[match(patients,rownames(cna)),]

    # merge #
    mut <- as.data.frame(cbind(mut,cna))
    rownames(mut) <- patients
  }

  if(spe.plat){
    v=strsplit(patients, "-IM")
    v=unlist(lapply(1:length(v), function(x)v[[x]][2]))

    # remove 410 platform patients #
    missing <- setdiff(paste0(g.impact$g468,c("",".fus",".Del",".Amp")),
                       paste0(g.impact$g410,c("",".fus",".Del",".Amp")))
    if(sum(v == "5") > 0 && sum(missing %in% colnames(mut)) > 0)
      mut[which(v == "5"), na.omit(match(missing, colnames(mut)))] <- NA

    # remove 341 platform patients #
    missing <- setdiff(paste0(g.impact$g468,c("",".fus",".Del",".Amp")),
                       paste0(g.impact$g341,c("",".fus",".Del",".Amp")))
    if(sum(v == "3") > 0 && sum(missing %in% colnames(mut)) > 0)
      mut[which(v == "3"), na.omit(match(missing, colnames(mut)))] <- NA
  }

  return(list("mut"=mut,"no.mut.patients"=rownames(mut)[missing.mut]))
}
