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
#' Default is NULL.
#' @param cna.binary A boolean argument specifying if the cna events should be enforced as binary. In which case separate columns for
#' amplifications and deletions will be created.
#' @param cna.relax for cna data only enables to count both gains and shallow deletions as amplifications and deletions respectively.
#' @param spe.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is TRUE.
#' @param set.plat character argument specifying which IMPACT platform the data should be reduced to if spe.plat is set to TRUE.
#'  Options are "341", "410" and "468". Default is NULL.
#' @param rm.empty boolean specifying if columns with no events founds should be removed. Default is TRUE.
#' @param pathway boolean specifying if pathway annotation should be applied. If TRUE, the function will return a supplementary binary
#' dataframe with columns being each pathway and each row being a sample. Default is FALSE.
#' @param col.names character vector of the necessary columns to be used. By default: col.names = c(Tumor_Sample_Barcode = NULL,
#'  Hugo_Symbol = NULL, Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL)
#' @param oncokb boolean specfiying if maf file should be oncokb annotated. Default is FALSE.
#' @param keep_onco A character vector specifying which oncoKB annotated variants to keep. Options are
#'  'Oncogenic', 'Likely Oncogenic', 'Predicted Oncogenic', 'Likely Neutral' and 'Inconclusive'. By default
#'  'Oncogenic', 'Likely Oncogenic' and 'Predicted Oncogenic' variants will be kept (recommended).
#' @param token the token affiliated to your oncoKB account.
#' @param ... Further arguments passed to the oncokb() function such a token
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
                   fusion = NULL,cna = NULL,cna.binary = TRUE,cna.relax = FALSE, spe.plat = TRUE,
                   set.plat = NULL,rm.empty = TRUE, pathway = FALSE,
                   col.names = c(Tumor_Sample_Barcode = NULL, Hugo_Symbol = NULL,
                                 Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL),
                   oncokb = FALSE, keep_onco = c("Oncogenic","Likely Oncogenic","Predicted Oncogenic"), token = '',...){

  impact_gene_info <- gnomeR::impact_gene_info

  if(is.null(maf) && is.null(fusion) && is.null(cna)) stop("You must provided one of the three following files: MAF, fusion or CNA.")
  # reformat columns #
  if(!is.null(maf)) {
    if("api" %in% class(maf)) is.api = TRUE
    maf <- as_tibble(maf) %>%
      rename(col.names)
  }
  if(!is.null(fusion)) fusion <- fusion %>%
      rename(col.names)

  # check oncokb API #
  if(oncokb)
    if(token == '')
      stop("If you want to OncoKB annotate your data you are required to have an API token.
           See https://www.oncokb.org/.")

  ## if data from API need to split mutations and fusions ##
  if(!is.null(maf)){
    if(is.na(match("Variant_Classification",colnames(maf))))
      stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")

    if(!is.null(maf) && is.null(fusion) &&
       nrow(as_tibble(maf) %>%
            filter(.data$Variant_Classification == "Fusion")) > 0){
      fusion <- as_tibble(maf) %>%
        filter(.data$Variant_Classification == "Fusion")
      if(is.api)
        fusion <- fusion %>%
          mutate(Fusion = gsub("fusion","",.data$proteinChange))

      maf <- as_tibble(maf) %>%
        filter(.data$Variant_Classification != "Fusion")
      warning("Fusions were found in the maf file, they were removed and a fusion file was created.")
    }
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

    # filter/define patients #
    if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
    else patients <- as.character(unique(maf$Tumor_Sample_Barcode))
    if(oncokb)
      maf <- oncokb(maf = maf, fusion = NULL, cna = NULL, token = token,...)$maf_oncokb %>%
        filter(.data$oncogenic %in% keep_onco)
    # set maf to maf class #
    maf <- structure(maf,class = c("data.frame","maf"))
    # getting mutation binary matrix #
    mut <- createbin(obj = maf, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)

  }

  # fusions #
  if(!is.null(fusion)){

    fusion <- as.data.frame(fusion)
    if(oncokb)
      fusion <- oncokb(maf = NULL, fusion = fusion, cna = NULL, token = token,...)$fusion_oncokb %>%
        filter(.data$oncogenic %in% keep_onco)
    fusion <- structure(fusion,class = c("data.frame","fusion"))
    # filter/define patients #
    if(is.null(patients)) patients <- as.character(unique(fusion$Tumor_Sample_Barcode))
    fusion <- createbin(obj = fusion, patients = patients, mut.type = mut.type, cna.binary = cna.binary,
                        SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,fusion))
      rownames(mut) <- patients}
    else mut <- fusion
  }

  # cna #
  if(!is.null(cna)){
    if("api" %in% class(cna)){

      ## oncokb for API ##
      if(oncokb){

        temp.cna <- as.data.frame(matrix(0L,
                                         nrow = length(unique(cna$Hugo_Symbol)),
                                         ncol = length(unique(cna$sampleId))+1))
        # rownames(temp.cna) <- unique(cna$SAMPLE_ID)
        temp.cna[,1] <- unique(cna$Hugo_Symbol)
        colnames(temp.cna) <- c("Hugo_Symbol",unique(cna$sampleId))

        for(i in colnames(temp.cna)[-1]){
          temp <- as_tibble(cna) %>%
            filter(.data$sampleId %in% i) %>%
            select(.data$sampleId, .data$Hugo_Symbol, .data$alteration)
          if(nrow(temp)>0){
            temp.cna[match(temp$Hugo_Symbol, temp.cna[,1]),match(i, colnames(temp.cna))] <- temp$alteration
          }
        }
        temp.cna[temp.cna == "Amplification"] <- 2
        temp.cna[temp.cna == "Deletion"] <- -2

        cna <- temp.cna
        temp.cna <- NULL

        cna <- oncokb(maf = NULL, fusion = NULL, cna = cna, token = token,...)$cna_oncokb %>%
          filter(.data$oncogenic %in% keep_onco) #%>%
        # dplyr::mutate(SAMPLE_ID = gsub("\\.","-",SAMPLE_ID))

        temp.cna <- as.data.frame(matrix(0L,
                                         nrow = length(unique(cna$HUGO_SYMBOL)),
                                         ncol = length(unique(cna$SAMPLE_ID))+1))
        # rownames(temp.cna) <- unique(cna$SAMPLE_ID)
        temp.cna[,1] <- unique(cna$HUGO_SYMBOL)
        colnames(temp.cna) <- c("Hugo_Symbol",unique(cna$SAMPLE_ID))

        for(i in colnames(temp.cna)[-1]){
          temp <- cna %>%
            filter(.data$SAMPLE_ID %in% i) %>%
            select(.data$SAMPLE_ID, .data$HUGO_SYMBOL, .data$ALTERATION)
          if(nrow(temp)>0){
            temp.cna[match(temp$HUGO_SYMBOL, temp.cna[,1]),match(i, colnames(temp.cna))] <- temp$ALTERATION
          }
        }
        temp.cna[temp.cna == "Amplification"] <- 2
        temp.cna[temp.cna == "Deletion"] <- -2

        cna <- temp.cna
        temp.cna <- NULL

        cna <- structure(cna,class = c("data.frame","cna"))
        if(is.null(patients)) patients <- gsub("\\.","-",as.character(colnames(cna)))[-1]
        cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                         SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)

      }

      else{
        if(is.null(patients)) patients <- unique(cna$sampleId)
        cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                         SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
      }
    }

    else{
      cna <- as.data.frame(cna)
      if(oncokb){
        cna <- oncokb(maf = NULL, fusion = NULL, cna = cna, token = token,...)$cna_oncokb %>%
          filter(.data$oncogenic %in% c("Oncogenic","Likely Oncogenic")) #%>%
          # dplyr::mutate(SAMPLE_ID = gsub("\\.","-",SAMPLE_ID))

        temp.cna <- as.data.frame(matrix(0L,
                                         nrow = length(unique(cna$HUGO_SYMBOL)),
                                         ncol = length(unique(cna$SAMPLE_ID))+1))
        # rownames(temp.cna) <- unique(cna$SAMPLE_ID)
        temp.cna[,1] <- unique(as.character(cna$HUGO_SYMBOL))
        colnames(temp.cna) <- c("Hugo_Symbol",unique(as.character(cna$SAMPLE_ID)))

        for(i in colnames(temp.cna)[-1]){
          temp <- cna %>%
            filter(.data$SAMPLE_ID %in% i) %>%
            select(.data$SAMPLE_ID, .data$HUGO_SYMBOL, .data$ALTERATION)
          if(nrow(temp)>0){
            temp.cna[match(temp$HUGO_SYMBOL, temp.cna[,1]),match(i, colnames(temp.cna))] <- as.character(temp$ALTERATION)
          }
        }
        temp.cna[temp.cna == "Amplification"] <- 2
        temp.cna[temp.cna == "Deletion"] <- -2

        cna <- temp.cna
        temp.cna <- NULL
      }

      cna <- structure(cna,class = c("data.frame","cna"))
      if(is.null(patients)) patients <- gsub("\\.","-",as.character(colnames(cna)))[-1]
      # else{
      #   colnames(cna) <- gsub("\\.","-",colnames(cna))
      # }
      cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                       SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    }
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

  # if(rm.empty && length(which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0))) mut <- mut[,which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0)]
  if(rm.empty && length(which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)))
    mut <- mut[,which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)]

  # create pathway levels alterations table #
  if(pathway){
    pathway_dat <- as.data.frame(do.call('cbind',lapply(unique(impact_gene_info$pathway[!is.na(impact_gene_info$pathway)]),function(x){
      genes <- as.character(unlist(impact_gene_info %>%
                                     filter(.data$pathway == x) %>% select(.data$hugo_symbol)))
      as.numeric(apply(mut %>% select(starts_with(genes)),1,function(y){
        ifelse(sum(abs(as.numeric(as.character(y))),na.rm = T)>0, 1,0)
      }))
    })))
    colnames(pathway_dat) <- unique(impact_gene_info$pathway[!is.na(impact_gene_info$pathway)])
    rownames(pathway_dat) <- rownames(mut)
    return(list(mut = mut, pathway_dat = pathway_dat))
  }

  return(mut)
}


##############################################
###### CREATE BINARIES FOR DIFF CLASSES ######
##############################################


createbin <- function(obj, patients, mut.type, cna.binary, SNP.only,include.silent, cna.relax, spe.plat){
  UseMethod("createbin")
}

createbin.default <- function(obj) {
  cat("The data did not match any known data type. Please review it and make sure it is correctly specified.")
}


##############################################
############# MUTATION MATRIX ################
##############################################

createbin.maf <- function(obj, patients, mut.type, cna.binary, SNP.only, include.silent, cna.relax, spe.plat){
  maf <- as_tibble(obj)
  maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
  if (sum(grepl("KMT2D", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        .data$Hugo_Symbol == "KMT2D" ~ "MLL2",
        TRUE ~ .data$Hugo_Symbol
      ))

    warning("KMT2D has been recoded to MLL2")
  }

  if (sum(grepl("KMT2C", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        .data$Hugo_Symbol == "KMT2C" ~ "MLL3",
        TRUE ~ .data$Hugo_Symbol
      ))

    warning("KMT2C has been recoded to MLL3")
  }

  if (sum(grepl("MYCL", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        .data$Hugo_Symbol == "MYCL" ~ "MYCL1",
        TRUE ~ .data$Hugo_Symbol
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

  maf <- as_tibble(maf) %>%
    filter(.data$Variant_Classification != Variant.filt,
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

createbin.fusion <- function(obj, patients, mut.type,cna.binary, SNP.only,include.silent, cna.relax, spe.plat){
  fusion <- as_tibble(obj)
  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a gene name column. (Hugo_Symbol)")


  fusion <- as_tibble(fusion) %>%
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

createbin.cna <- function(obj, patients, mut.type,cna.binary, SNP.only,include.silent, cna.relax, spe.plat){
  cna <- obj
  rownames(cna) <- cna[,1]
  cna <- cna[,-1]
  cna <- as.data.frame(t(cna))
  rownames(cna) <- gsub("\\.","-",rownames(cna))
  cna <- cna[rownames(cna) %in% patients,]

  if(cna.binary){
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
    if(length(which(is.na(match(patients,rownames(cna))))) > 0){
      missing <- patients[which(is.na(match(patients,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
      cna <- cna[match(patients,rownames(cna)),]
    }
  }
  if(!cna.binary){
    # add missing
    if(length(which(is.na(match(patients,rownames(cna))))) > 0){
      missing <- patients[which(is.na(match(patients,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
      cna <- cna[match(patients,rownames(cna)),]
    }

    cna <- cna %>%
      mutate_all(~ factor(as.numeric(as.character(.)),
                          levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(.)))]))
    colnames(cna) <- paste0(colnames(cna),".cna")
  }

  return(cna)
}

### cna from API ###
createbin.api <- function(obj, patients, mut.type,cna.binary, SNP.only,include.silent, cna.relax, spe.plat){
  cna <- as.data.frame(obj)

  # recreate orginal format #
  temp <- as.data.frame(matrix(0L,ncol = length(patients)+1, nrow = length(unique(cna$Hugo_Symbol))))
  colnames(temp) <- c("Hugo_Symbol",patients)
  temp[,1] <- unique(cna$Hugo_Symbol)
  for(i in patients){
    temp[match(as.character(unlist(cna %>% filter(.data$sampleId %in% i) %>% select(.data$Hugo_Symbol))),temp[,1]),
         match(i, colnames(temp))] <- as.numeric(unlist(cna %>% filter(.data$sampleId %in% i) %>% select(.data$alteration)))
  }

  cna <- temp
  rownames(cna) <- cna[,1]
  cna <- cna[,-1]
  cna <- as.data.frame(t(cna))
  cna <- cna[rownames(cna) %in% patients,]

  if(cna.binary){
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
    if(length(which(is.na(match(patients,rownames(cna))))) > 0){
      missing <- patients[which(is.na(match(patients,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
      cna <- cna[match(patients,rownames(cna)),]
    }
  }
  if(!cna.binary){
    # add missing
    if(length(which(is.na(match(patients,rownames(cna))))) > 0){
      missing <- patients[which(is.na(match(patients,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
      cna <- cna[match(patients,rownames(cna)),]
    }

    cna <- cna %>%
      mutate_all(~ factor(as.numeric(as.character(.)),
                          levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(.)))]))
    colnames(cna) <- paste0(colnames(cna),".cna")
  }

  return(cna)
}
