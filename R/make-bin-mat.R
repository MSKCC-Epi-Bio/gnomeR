#' binmat
#' Enables creation of a binary matrix from a maf file with
#' a predefined list of patients (rows are patients and columns are genes)
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
#' Default is NULL. Note that this file must have patients as columns and genes as rows. binmat expects a matrix with
#' values between -2 and 2. Please do not use any other format. Other functions in the package are available to deal with more detailed
#' CNA data.
#' @param cna.relax for cna data only enables to count both gains and shallow deletions as amplifications and deletions respectively.
#' @param spe.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is TRUE.
#' @param set.plat character argument specifying which IMPACT platform the data should be reduced to if spe.plat is set to TRUE.
#'  Options are "341", "410" and "468". Default is NULL.
#' @param rm.empty boolean specifying if columns with no events founds should be removed. Default is TRUE.
#' @param col.names character vector of the necessary columns to be used. By default: col.names = c(Tumor_Sample_Barcode = NULL,
#'  Hugo_Symbol = NULL, Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL)
#' @return mut : a binary matrix of mutation data
#' @export
#' @examples library(gnomeR)
#' # mut.only <- binmat(maf = mut)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,maf = mut,
#' mut.type = "SOMATIC",SNP.only = FALSE,
#' include.silent = FALSE, spe.plat = FALSE)
#' bin.mut <- binmat(patients = patients,maf = mut,
#' mut.type = "SOMATIC",SNP.only = FALSE,
#' include.silent = FALSE,
#' cna.relax = TRUE, spe.plat = FALSE,
#'  set.plat = "410", rm.empty = FALSE)
#' @import dplyr
#' @import dtplyr
#' @import stringr


###############################################
###### MAIN FUNCTION GROUPING EVERYTHING ######
###############################################

binmat <- function(patients=NULL, maf = NULL, mut.type = "SOMATIC",SNP.only = FALSE,include.silent = FALSE,
                   fusion = NULL,cna = NULL,cna.relax = FALSE, spe.plat = TRUE, set.plat = NULL,rm.empty = TRUE,
                   col.names = c(Tumor_Sample_Barcode = NULL, Hugo_Symbol = NULL,
                                 Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL)){

  if(is.null(maf) && is.null(fusion) && is.null(cna)) stop("You must provided one of the three following files: MAF, fusion or CNA.")
  # reformat columns #
  if(!is.null(maf)) maf <- maf %>%
      rename(col.names)
  if(!is.null(fusion)) fusion <- fusion %>%
      rename(col.names)

  ## if data from API need to split mutations and fusions ##
  if(!is.null(maf) && is.null(fusion) &&
     nrow(maf %>%
          filter(.data$Variant_Classification == "Fusion")) > 0){
    fusion <- maf %>%
      filter(.data$Variant_Classification == "Fusion")
    maf <- maf %>%
      filter(.data$Variant_Classification != "Fusion")
    warning("Fusions were found in the maf file, they were removed and a fusion file was created.")
  }


  mut <- NULL

  if(!is.null(maf)){

    # quick data checks #
    if(is.na(match("Tumor_Sample_Barcode",colnames(maf))))
      stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
    if(is.na(match("Hugo_Symbol",colnames(maf))))
      stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
    if(is.na(match("Variant_Classification",colnames(maf))))
      stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
    if(is.na(match("Mutation_Status",colnames(maf)))){
      warning("The MAF file inputted is missing a mutation status column (Mutation_Status). It will be assumed that
            all variants are of the same type (SOMATIC/GERMLINE).")
      maf$Mutation_Status <- rep("SOMATIC",nrow(maf))
    }
    if(is.na(match("Variant_Type",colnames(maf)))){
      warning("The MAF file inputted is missing a mutation status column (Variant_Type). It will be assumed that
            all variants are of the same type (SNPs).")
      maf$Variant_Type <- rep("SNPs",nrow(maf))
    }

    # set maf to maf class #
    maf <- structure(maf,class = c("data.frame","maf"))
    # filter/define patients #
    if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
    else patients <- as.character(unique(maf$Tumor_Sample_Barcode))
    # getting mutation binary matrix #
    mut <- createbin(obj = maf, patients = patients, mut.type = mut.type,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)

  }

  # fusions #
  if(!is.null(fusion)){

    fusion <- as.data.frame(fusion)
    fusion <- structure(fusion,class = c("data.frame","fusion"))
    # filter/define patients #
    if(is.null(patients)) patients <- as.character(unique(fusion$Tumor_Sample_Barcode))
    fusion <- createbin(obj = fusion, patients = patients, mut.type = mut.type,
                        SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,fusion))
      rownames(mut) <- patients}
    else mut <- fusion
  }

  # cna #
  if(!is.null(cna)){
    cna <- as.data.frame(cna)
    cna <- structure(cna,class = c("data.frame","cna"))
    if(is.null(patients)) patients <- as.character(colnames(cna))
    cna <- createbin(obj = cna, patients = patients, mut.type = mut.type,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,cna))
      rownames(mut) <- patients}
    else mut <- cna
  }

  # specific platform for IMPACT #
  if(spe.plat){

    v=strsplit(patients, "-IM|-IH")
    if(!all(lapply(v, length) == 2)){
      warning("All patients were not sequenced on the IMPACT platform or some were mispecified. '-IM' or '-IH' requiered in sample ID.
              The spe.plat argument has been overwritten to FALSE.")
      spe.plat = F
    }
    v=unlist(lapply(1:length(v), function(x)v[[x]][2]))
    if(length(unique(v)) == 1){
      warning("All samples were sequenced on the same platform.
              The spe.plat argument has been overwritten to FALSE.")
      spe.plat = F
    }
    if(spe.plat){
      g.impact <- g.impact
      # remove 410 platform patients #
      missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp")),
                         c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp")))
      if(sum(v == "5") > 0 && sum(missing %in% colnames(mut)) > 0)
        mut[which(v == "5"), stats::na.omit(match(missing, colnames(mut)))] <- NA

      # remove 341 platform patients #
      missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp")),
                         c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp")))
      if(sum(v == "3") > 0 && sum(missing %in% colnames(mut)) > 0)
        mut[which(v == "3"), stats::na.omit(match(missing, colnames(mut)))] <- NA

    }
  }

  if(!is.null(set.plat)){
    if(set.plat == "341"){
      keep <- c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp"))
      mut <- mut[,colnames(mut) %in% keep]
    }
    if(set.plat == "410"){
      keep <- c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp"))
      mut <- mut[,colnames(mut) %in% keep]
    }
  }

  if(rm.empty && length(which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0))) mut <- mut[,which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0)]
  return(mut)
}


##############################################
###### CREATE BINARIES FOR DIFF CLASSES ######
##############################################


createbin <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  UseMethod("createbin")
}

createbin.default <- function(obj) {
  cat("The data did not match any known data type. Please review it and make sure it is correctly specified.")
}


##############################################
############# MUTATION MATRIX ################
##############################################

createbin.maf <- function(obj, patients, mut.type, SNP.only, include.silent, cna.relax, spe.plat){
  maf <- obj
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

  if (sum(grepl("MYCL", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "MYCL" ~ "MYCL1",
        TRUE ~ Hugo_Symbol
      ))

    warning("MYCL has been recoded to MYCL1")
  }

  # # clean gen dat #
  if(SNP.only) SNP.filt = "SNP"
  else SNP.filt = unique(maf$Variant_Type)

  if(!include.silent) Variant.filt = "Silent"
  else Variant.filt = ""

  if(tolower(mut.type) == "all") Mut.filt = unique(maf$Mutation_Status)
  else Mut.filt = mut.type

  maf <- maf %>% filter(.data$Variant_Classification != Variant.filt,
                        .data$Variant_Type %in% SNP.filt,
                        tolower(.data$Mutation_Status) %in% tolower(Mut.filt))


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
    warning(paste0("Some patients did not have any mutations found in the MAF file.", paste0(rownames(mut)[missing.mut], collapse = ",")))

  return(mut)
}


###########################################
############# FUSION MATRIX ###############
###########################################

createbin.fusion <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  fusion <- obj
  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a gene name column. (Hugo_Symbol)")


  fusion <- fusion %>%
    filter(.data$Tumor_Sample_Barcode %in% patients)

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
  return(fusion.out)
}


################################################
############# COPY NUMBER MATRIX ###############
################################################

createbin.cna <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  cna <- obj
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

  return(cna)
}
