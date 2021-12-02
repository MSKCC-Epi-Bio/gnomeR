#' binmat
#' \%lifecycle{stable}
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
#' @param cna.relax By default this argument is set to FALSE, where only deep deletions (-2) and amplifications (2) will be annotated as events.
#'  When set to FTRUE all deletions (-1 shallow and -2 deep) are counted as an event same for all gains (1 gain, 2 amplification) as an event.
#' @param specify.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is TRUE.
#' @param set.plat character argument specifying which IMPACT platform the data should be reduced to if specify.plat is set to TRUE.
#'  Options are "341", "410" and "468. Default is NULL.
#' @param rm.empty boolean specifying if columns with no events founds should be removed. Default is TRUE.
#' @param pathway boolean specifying if pathway annotation should be applied. If TRUE, the function will return a supplementary binary
#' dataframe with columns being each pathway and each row being a sample. Default is FALSE.
#' @param recode.aliases bolean specifying if automated gene name alias matching should be done. Default is TRUE. When TRUE
#' the function will check for genes that may have more than 1 name in your data using the aliases im gnomeR::impact_gene_info alias column
#' @param col.names character vector of the necessary columns to be used. By default: col.names = c(Tumor_Sample_Barcode = NULL,
#'  Hugo_Symbol = NULL, Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL)
#' @param oncokb boolean specfiying if maf file should be oncokb annotated. Default is FALSE.
#' @param keep_onco A character vector specifying which oncoKB annotated variants to keep. Options are
#'  'Oncogenic', 'Likely Oncogenic', 'Predicted Oncogenic', 'Likely Neutral' and 'Inconclusive'. By default
#'  'Oncogenic', 'Likely Oncogenic' and 'Predicted Oncogenic' variants will be kept (recommended).
#' @param token the token affiliated to your oncoKB account.
#' @param sample_panels A dataframe containing the sample ids corresponding to the patients argument in the first
#' column and their corresponding genomic panel. At the moment gnomeR can handle panels from the following centers:
#' MSKCC, DFCI, VICC and UHN. See objects impact_gene_info and panel_names for more information.
#' @param ... Further arguments passed to the oncokb() function such a token
#' @return mut : a binary matrix of mutation data
#' @export
#' @examples library(gnomeR)
#' # mut.only <- binmat(maf = mut)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,maf = mut,
#' mut.type = "SOMATIC",SNP.only = FALSE,
#' include.silent = FALSE, specify.plat = FALSE)
#' bin.mut <- binmat(patients = patients,maf = mut,
#' mut.type = "SOMATIC",SNP.only = FALSE,
#' include.silent = FALSE,
#' cna.relax = TRUE, specify.plat = FALSE,
#'  set.plat = "410", rm.empty = FALSE)
#' @import dplyr
#' @import dtplyr
#' @import stringr


###############################################
###### MAIN FUNCTION GROUPING EVERYTHING ######
###############################################

binmat <- function(patients=NULL, maf = NULL, mut.type = "SOMATIC",SNP.only = FALSE,include.silent = FALSE,
                   fusion = NULL,cna = NULL,cna.binary = TRUE,cna.relax = FALSE, specify.plat = TRUE,
                   set.plat = NULL,rm.empty = TRUE, pathway = FALSE,
                   recode.aliases = TRUE,
                   col.names = c(Tumor_Sample_Barcode = NULL, Hugo_Symbol = NULL,
                                 Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL),
                   oncokb = FALSE, keep_onco = c("Oncogenic","Likely Oncogenic","Predicted Oncogenic"),
                   token = '', sample_panels = NULL,...){

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

  # make patients/samples list #
  if(is.null(patients)){
    patients <- c()
    if(!is.null(maf))
      patients <- unique(c(patients, as.character(unique(maf$Tumor_Sample_Barcode))))
    if(!is.null(fusion))
      patients <- unique(c(patients, as.character(unique(fusion$Tumor_Sample_Barcode))))
    if(!is.null(cna))
      if("api" %in% class(cna))
        patients <- unique(c(patients, as.character(unique(cna$sampleId))))
      else if(length(grep("oncogenic", colnames(cna), ignore.case = TRUE)) == 0)
        patients <- unique(c(patients, as.character(unique(gsub("\\.","-",colnames(cna)[-1])))))
      else if(length(grep("oncogenic", colnames(cna), ignore.case = TRUE)) > 0)
        patients <- unique(c(patients, as.character(unique(cna$SAMPLE_ID))))
  }
  else
    patients <- as.character(patients)

  mut <- NULL

  if(!is.null(maf)){

    # quick data checks #
    maf <- check_maf_input(maf, recode.aliases= recode.aliases)
    # if(is.na(match("Tumor_Sample_Barcode",colnames(maf))))
    #   stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
    # if(is.na(match("Hugo_Symbol",colnames(maf))))
    #   stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
    # if(is.na(match("Variant_Classification",colnames(maf))))
    #   stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
    # if(is.na(match("Mutation_Status",colnames(maf)))){
    #   warning("The MAF file inputted is missing a mutation status column (Mutation_Status). It will be assumed that
    #         all variants are of the same type (SOMATIC/GERMLINE).")
    #   maf$Mutation_Status <- rep("SOMATIC",nrow(maf))
    # }
    # if(is.na(match("Variant_Type",colnames(maf)))){
    #   warning("The MAF file inputted is missing a mutation status column (Variant_Type). It will be assumed that
    #         all variants are of the same type (SNPs).")
    #   maf$Variant_Type <- rep("SNPs",nrow(maf))
    # }

    # filter/define patients #
    # if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
    # else patients <- as.character(unique(maf$Tumor_Sample_Barcode))
    if(oncokb)
      maf <- oncokb(maf = maf, fusion = NULL, cna = NULL, token = token,...)$maf_oncokb %>%
        filter(.data$ONCOGENIC %in% keep_onco)
    # set maf to maf class #
    maf <- structure(maf,class = c("data.frame","maf"))
    # getting mutation binary matrix #
    mut <- createbin(obj = maf, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, specify.plat = specify.plat, recode.aliases = recode.aliases)

  }

  # fusions #
  if(!is.null(fusion)){

    fusion <- as.data.frame(fusion)
    if(oncokb)
      fusion <- oncokb(maf = NULL, fusion = fusion, cna = NULL, token = token,...)$fusion_oncokb %>%
        filter(.data$ONCOGENIC %in% keep_onco)
    fusion <- structure(fusion,class = c("data.frame","fusion"))
    # filter/define patients #
    # if(is.null(patients)) patients <- as.character(unique(fusion$Tumor_Sample_Barcode))
    fusion <- createbin(obj = fusion, patients = patients, mut.type = mut.type, cna.binary = cna.binary,
                        SNP.only = SNP.only, include.silent = include.silent, specify.plat = specify.plat, recode.aliases = recode.aliases)
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
        cna <- oncokb(maf = NULL, fusion = NULL, cna = cna, token = token,...)$cna_oncokb %>%
          filter(.data$ONCOGENIC %in% keep_onco) #%>%
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
        # if(is.null(patients)) patients <- gsub("\\.","-",as.character(colnames(cna)))[-1]
        cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                         SNP.only = SNP.only, include.silent = include.silent, specify.plat = specify.plat, recode.aliases = recode.aliases)

      }

      else{
        # if(is.null(patients)) patients <- unique(cna$sampleId)
        cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                         SNP.only = SNP.only, include.silent = include.silent, specify.plat = specify.plat,
                         recode.aliases = recode.aliases)
      }
    }

    else{
      cna <- as.data.frame(cna)

      if(length(grep("oncogenic",colnames(cna), ignore.case = TRUE)) > 0){
        if(oncokb){
          warning("Your data is already oncokb annotated. 'oncokb' argument was overwriten to FALSE.")
          oncokb <- FALSE
        }

        cna <- cna %>%
          filter(.data$ONCOGENIC %in% keep_onco) #%>%
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


      if(oncokb){
        cna <- oncokb(maf = NULL, fusion = NULL, cna = cna, token = token,...)$cna_oncokb %>%
          filter(.data$ONCOGENIC %in% keep_onco) #%>%
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
      # if(is.null(patients)) patients <- gsub("\\.","-",as.character(colnames(cna)))[-1]
      # else{
      #   colnames(cna) <- gsub("\\.","-",colnames(cna))
      # }
      cna <- createbin(obj = cna, patients = patients, mut.type = mut.type, cna.binary = cna.binary,cna.relax = cna.relax,
                       SNP.only = SNP.only, include.silent = include.silent, specify.plat = specify.plat, recode.aliases = recode.aliases)
    }
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,cna))
      rownames(mut) <- patients}
    else mut <- cna
  }

  # fix colnames #
  colnames(mut) <- gsub("-",".",colnames(mut))


  ## ANNOTATE SAMPLES ##
  if(specify.plat){
    if(is.null(sample_panels)){
      warning("When panels argument is NULL only IMPACT samples can be annotated for missingness.
          If you are using GENIE data, we recommend that you provide this argument.")

      ## ANNOTATE ONLY IMPACT SAMPLES ##
      inds <- grep("-IM|-IH",rownames(mut))
      if(length(inds) > 0){
        mut_to_annotate <- mut[inds,]
        panels_to_use <- rownames(mut)[inds]

        panels_to_use <- unique(sapply(strsplit(panels_to_use,"-IM|-IH"), FUN = function(x){x[2]}))

        mut_annotated <- as.data.frame(
          do.call('rbind',
                  lapply(panels_to_use, function(p){

                    if(p == "3"){
                      p_name <- "MSK_341"
                      samples_to_annotate <-
                        as.character(rownames(mut_to_annotate)[grep("-IM3|-IH3",rownames(mut_to_annotate))])
                    }
                    if(p == "5"){
                      p_name <- "MSK_410"
                      samples_to_annotate <-
                        as.character(rownames(mut_to_annotate)[grep("-IM5|-IH5",rownames(mut_to_annotate))])
                    }
                    # need to add IM7 specific panel #
                    if(p %in% c("6","7")){
                      p_name <- "MSK_468"
                      samples_to_annotate <-
                        as.character(rownames(mut_to_annotate)[grep("-IM6|-IH6",rownames(mut_to_annotate))])
                    }
                    # need to change this to all genes that are not included! #
                    gene_to_NA <- setdiff(gsub(".fus|.Del|.Amp|.cna","",colnames(mut)) ,
                                          as.character(unlist(genie_gene_info %>%
                                                                select(hugo_symbol,gsub("-","_",p_name)) %>%
                                                                filter(!!as.symbol(p_name) == "included") %>%
                                                                select(hugo_symbol))))

                    missing <- c(gene_to_NA, paste0(gene_to_NA,".fus"),paste0(gene_to_NA,".Del"),
                                 paste0(gene_to_NA,".Amp"),paste0(gene_to_NA,".cna"))

                    mut_sub <- mut[samples_to_annotate,]
                    mut_sub[,na.omit(match(missing,colnames(mut_sub)))] <- NA
                    return(mut_sub)
                  })
          )
        )

        mut <- rbind(mut[-inds,], mut_annotated)
      }
    }
    else{
      panels_to_use <- unique(as.character(sample_panels[,2]))
      mut_annotated <- as.data.frame(
        do.call('rbind',
                lapply(panels_to_use, function(p){
                  p_name <- panel_names$Panel[match(p, panel_names$Sequence.Assay.ID)]
                  samples_to_annotate <- as.character(sample_panels[sample_panels[,2] == p,1])
                  # need to change this to all genes that are not included! #
                  gene_to_NA <- setdiff(gsub(".fus|.Del|.Amp|.cna","",colnames(mut)) ,
                                        as.character(unlist(genie_gene_info %>%
                                                              select(hugo_symbol,gsub("-","_",p_name)) %>%
                                                              filter(!!as.symbol(p_name) == "included") %>%
                                                              select(hugo_symbol))))

                  missing <- c(gene_to_NA, paste0(gene_to_NA,".fus"),paste0(gene_to_NA,".Del"),
                               paste0(gene_to_NA,".Amp"),paste0(gene_to_NA,".cna"))

                  mut_sub <- mut[samples_to_annotate,]
                  mut_sub[,na.omit(match(missing,colnames(mut_sub)))] <- NA
                  return(mut_sub)
                })
        )
      )
      mut <- mut_annotated
    }
  }


  # define g.impact from impact_gene_info
  g.impact <- list()
  g.impact$g341 <- genie_gene_info$hugo_symbol[genie_gene_info$MSK_341 == "included"]
  g.impact$g410 <- genie_gene_info$hugo_symbol[genie_gene_info$MSK_410 == "included"]
  g.impact$g468 <- genie_gene_info$hugo_symbol[genie_gene_info$MSK_468 == "included"]


  if(!is.null(set.plat)){
    if(set.plat == "341"){
      keep <- c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp"),paste0(g.impact$g341,".cna"))
      mut <- mut[,colnames(mut) %in% keep]
    }
    if(set.plat == "410"){
      keep <- c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp"),paste0(g.impact$g410,".cna"))
      mut <- mut[,colnames(mut) %in% keep]
    }
    if(set.plat == "468"){
      keep <- c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp"),paste0(g.impact$g468,".cna"))
      mut <- mut[,colnames(mut) %in% keep]
    }
  }


  # if(rm.empty && length(which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0))) mut <- mut[,which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0)]
  if(rm.empty && length(which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)))
    mut <- mut[,which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)]




  # define g.impact from impact_gene_info
  # g.impact <- list()
  # g.impact$g341 <- impact_gene_info$hugo_symbol[impact_gene_info$platform_341 == "included"]
  # g.impact$g410 <- impact_gene_info$hugo_symbol[impact_gene_info$platform_410 == "included"]
  # g.impact$g468 <- impact_gene_info$hugo_symbol[impact_gene_info$platform_468 == "included"]
  #
  # # specific platform for IMPACT #
  # if(specify.plat){
  #
  #   v=strsplit(patients, "-IM|-IH")
  #   if(!all(lapply(v, length) == 2)){
  #     warning("All patients were not sequenced on the IMPACT platform or some were mispecified.
  #             Only samples endind in '-IM' or '-IH' will be annotated for specific IMPACT platforms.")
  #     # specify.plat = F
  #   }
  #   v=unlist(lapply(1:length(v), function(x)v[[x]][2]))
  #   if(length(unique(v[!is.na(v)])) == 1){
  #     warning("All samples were sequenced on the same platform.
  #             The specify.plat argument has been overwritten to FALSE.")
  #     specify.plat = F
  #   }
  #   if(specify.plat){
  #
  #     # remove 410 platform patients #
  #     missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp"),paste0(g.impact$g468,".cna")),
  #                        c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp"),paste0(g.impact$g410,".cna")))
  #     if(sum(v == "5", na.rm = T) > 0 && sum(missing %in% colnames(mut)) > 0)
  #       mut[which(v == "5"), stats::na.omit(match(missing, colnames(mut)))] <- NA
  #
  #     # remove 341 platform patients #
  #     missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp"),paste0(g.impact$g468,".cna")),
  #                        c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp"),paste0(g.impact$g341,".cna")))
  #     if(sum(v == "3", na.rm = T) > 0 && sum(missing %in% colnames(mut)) > 0)
  #       mut[which(v == "3"), stats::na.omit(match(missing, colnames(mut)))] <- NA
  #
  #   }
  # }
  #
  # if(!is.null(set.plat)){
  #   if(set.plat == "341"){
  #     keep <- c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp"),paste0(g.impact$g341,".cna"))
  #     mut <- mut[,colnames(mut) %in% keep]
  #   }
  #   if(set.plat == "410"){
  #     keep <- c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp"),paste0(g.impact$g410,".cna"))
  #     mut <- mut[,colnames(mut) %in% keep]
  #   }
  # }
  #
  # # if(rm.empty && length(which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0))) mut <- mut[,which(apply(mut,2,function(x){sum(x,na.rm=TRUE)})>0)]
  # if(rm.empty && length(which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)))
  #   mut <- mut[,which(apply(mut,2,function(x){length(unique(x[!is.na(x)]))})>1)]

  # create pathway levels alterations table #
  if(pathway){

    pathway_dat <- purrr::map2_dfc(pathways$pathway, pathways$genes, function(x, y) {

      genes <- y[y %in% names(mut)]
      res <- apply(mut %>% select(genes), 1,
                   function(y){
                     ifelse(sum(abs(as.numeric(as.character(y))),
                                na.rm = T) > 0, 1, 0)}) %>%
        as.numeric(.) %>%
        as.data.frame()

      names(res) <- x
      return(res)
    }
    )

    rownames(pathway_dat) <- rownames(mut)

    mut <- list(mut = mut, pathway_dat = pathway_dat)
  }
  # get possible genes that were not part of IMPACT #
  unique_genes <- unique(gsub(".Del|.Amp|.fus|.cna","",colnames(mut)))
  missing_genes <- unique_genes[which(!(unique_genes %in% genie_gene_info$hugo_symbol))]
  if(length(missing_genes) > 0 && specify.plat)
    warning(paste0("The following genes in the final matrix were not part of the official IMPACT or GENIE panel and
    thus couldn't be annotated for missing status.
                   To see a complete list of genes in IMPACT and GENIE please see 'genie_gene_info': ",
                   paste0(missing_genes,collapse = ", ")
    )
    )
  mut
}


##############################################
###### CREATE BINARIES FOR DIFF CLASSES ######
##############################################


createbin <- function(obj, patients, mut.type, cna.binary, SNP.only,include.silent, cna.relax, specify.plat, recode.aliases){
  UseMethod("createbin")
}

createbin.default <- function(obj) {
  cat("The data did not match any known data type. Please review it and make sure it is correctly specified.")
}


##############################################
############# MUTATION MATRIX ################
##############################################

createbin.maf <- function(obj, patients, mut.type, cna.binary, SNP.only, include.silent, cna.relax, specify.plat, recode.aliases = recode.aliases){
  maf <- as_tibble(obj)
  maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)

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
    warning(paste0("Some patients did not have any mutations found in the MAF file. (",
                   sum(missing.mut), ",",
                   round(sum(missing.mut)/nrow(mut)*100, digits = 2),"%): ",
                   paste0(rownames(mut)[missing.mut], collapse = ",")))

  return(mut)
}


###########################################
############# FUSION MATRIX ###############
###########################################

createbin.fusion <- function(obj, patients, mut.type,cna.binary,
                             SNP.only,include.silent, cna.relax,
                             specify.plat, recode.aliases){
  fusion <- as_tibble(obj)
  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a gene name column. (Hugo_Symbol)")
  fusion$Hugo_Symbol <- as.character(fusion$Hugo_Symbol)

  # get table of gene aliases
  alias_table <- tidyr::unnest(impact_gene_info, cols = alias) %>%
    select(hugo_symbol, alias)

  # recode aliases ---
  if(recode.aliases == TRUE) {

    fusion$Hugo_Symbol_Old <- fusion$Hugo_Symbol
    fusion$Hugo_Symbol <- purrr::map_chr(fusion$Hugo_Symbol, ~resolve_alias(.x,
                                                                            alias_table = alias_table))

    message <- fusion %>%
      dplyr::filter(Hugo_Symbol_Old != Hugo_Symbol) %>%
      dplyr::select(Hugo_Symbol_Old, Hugo_Symbol) %>%
      dplyr::distinct()

    if(nrow(message) > 0) {
      warning(paste0("FUSION DATA: To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your fusion data have been recoded. You can supress this with recode.aliases = FALSE \n \n",
                     purrr::map2(message$Hugo_Symbol_Old,
                                 message$Hugo_Symbol,
                                 ~paste0(.x, " recoded to ", .y, " \n"))))
    }
  }



  fusion <- as_tibble(fusion) %>%
    filter(.data$Tumor_Sample_Barcode %in% patients)

  #### out frame
  fusion.out <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(fusion$Hugo_Symbol))))
  colnames(fusion.out) <- unique(fusion$Hugo_Symbol)
  rownames(fusion.out) <- patients

  for(i in patients){
    genes <- fusion$Hugo_Symbol[fusion$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0)
      fusion.out[match(i,rownames(fusion.out)),
                 match(unique(as.character(genes)),colnames(fusion.out))] <- 1

  }
  colnames(fusion.out) <- paste0(colnames(fusion.out),".fus")
  return(fusion.out)
}


################################################
############# COPY NUMBER MATRIX ###############
################################################

createbin.cna <- function(obj, patients, mut.type,cna.binary,
                          SNP.only,include.silent, cna.relax,
                          specify.plat, recode.aliases){
  cna <- obj
  cna <- as.data.frame(tibble::as_tibble(cna))
  cna$Hugo_Symbol <- as.character(cna$Hugo_Symbol)

  # get table of gene aliases
  alias_table <- tidyr::unnest(impact_gene_info, cols = alias) %>%
    select(hugo_symbol, alias)

  # recode aliases ---
  if(recode.aliases == TRUE) {

    cna$Hugo_Symbol_Old <- cna$Hugo_Symbol
    cna$Hugo_Symbol <- purrr::map_chr(cna$Hugo_Symbol, ~resolve_alias(.x,
                                                                      alias_table = alias_table))

    message <- cna %>%
      dplyr::filter(Hugo_Symbol_Old != Hugo_Symbol) %>%
      dplyr::select(Hugo_Symbol_Old, Hugo_Symbol) %>%
      dplyr::distinct()

    if(nrow(message) > 0) {
      warning(paste0("CNA DATA: To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your CNA data have been recoded. You can supress this with recode.aliases = FALSE. \n \n",
                     purrr::map2(message$Hugo_Symbol_Old,
                                 message$Hugo_Symbol,
                                 ~paste0(.x, " recoded to ", .y, " \n"))))
    }
    cna <- cna %>% select(-one_of("Hugo_Symbol_Old"))
  }


  dups <- cna$Hugo_Symbol[duplicated(cna$Hugo_Symbol)]
  if(length(dups) > 0){
    for(i in dups){
      temp <- cna[which(cna$Hugo_Symbol == i),] # grep(i, cna$Hugo_Symbol,fixed = TRUE)
      temp2 <- as.character(unlist(apply(temp, 2, function(x){
        if(all(is.na(x)))
          out <- NA
        else if(anyNA(x))
          out <- x[!is.na(x)]
        else if(length(unique(x)) > 1)
          out <- x[which(x != 0)]
        else
          out <- x[1]
        return(out)
      })))
      temp2[-1] <- as.numeric(temp2[-1])
      cna <- rbind(cna[-which(cna$Hugo_Symbol == i),],
                   temp2)
    }
  }

  rownames(cna) <- cna$Hugo_Symbol
  cna <- cna[,-1]
  cna <- as.data.frame(t(cna))
  rownames(cna) <- gsub("\\.","-",rownames(cna))
  cna <- cna[rownames(cna) %in% patients,]
  samples_temp <- rownames(cna)

  if(cna.binary){
    temp <- do.call("cbind",apply(cna,2,function(x){
      if(cna.relax){
        yA <- ifelse(as.numeric(x)>=0.9,1,0)
        yD <- ifelse(as.numeric(x)<=-0.9,1,0)
      }
      if(!cna.relax){
        yA <- ifelse(as.numeric(x)==2,1,0)
        yD <- ifelse(as.numeric(x)<=-0.9,1,0) #==-2
      }
      out <- as.data.frame(cbind(yA,yD))
      colnames(out) <- c("Amp","Del")
      return(out)
    }))

    cna <- temp[,apply(temp,2,function(x){sum(x,na.rm=T) > 0})]
    rownames(cna) <- samples_temp

    # add missing
    if(length(which(is.na(match(patients,rownames(cna))))) > 0){
      missing <- patients[which(is.na(match(patients,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
    }
    cna <- cna[match(patients,rownames(cna)),]
    cna[is.na(cna)] <- 0
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
    cna <- cna[match(patients,rownames(cna)),]

    cna <- cna %>%
      mutate_all(~ as.numeric(gsub(" ","",as.character(.)))) %>%
      mutate_all(
        ~ case_when(
          . == 0 ~ "NEUTRAL",
          . %in% c(-1.5,-1) ~ "LOH",
          . == 1 ~ "GAIN",
          . == 2 ~ "AMPLIFICATION",
          . == -2 ~ "DELETION",
        )
      ) %>%
      mutate_all(
        ~factor(.,
                levels = c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION")[
                  which(c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION") %in% .)
                ])
      )

    # cna <- cna %>%
    #   mutate_all(~ factor(as.numeric(as.character(.)),
    #                       levels = c("0","-2","-1.5","2")[which(c(0,-2,-1.5,2) %in% as.numeric(as.character(.)))]))
    colnames(cna) <- paste0(colnames(cna),".cna")
    rownames(cna) <- patients

    cna[is.na(cna)] <- "NEUTRAL"
  }

  # cna[is.na(cna)] <- 0
  return(cna)
}

### cna from API ###
createbin.api <- function(obj, patients, mut.type,cna.binary,
                          SNP.only,include.silent, cna.relax, specify.plat, recode.aliases){
  cna <- as.data.frame(obj)

  cna <- as.data.frame(tibble::as_tibble(cna))
  cna$Hugo_Symbol <- as.character(cna$Hugo_Symbol)

  # get table of gene aliases
  alias_table <- tidyr::unnest(impact_gene_info, cols = alias) %>%
    select(hugo_symbol, alias)

  # recode aliases
  # cna$Hugo_Symbol_Old <- cna$Hugo_Symbol
  cna$Hugo_Symbol <- purrr::map_chr(cna$Hugo_Symbol, ~resolve_alias(.x,
                                                                    alias_table = alias_table))


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
  samples_temp <- rownames(cna)

  if(cna.binary){
    temp <- do.call("cbind",apply(cna,2,function(x){
      if(cna.relax){
        yA <- ifelse(x>=0.9,1,0)
        yD <- ifelse(as.numeric(x)<=-0.9,1,0) #==-2 # yD <- ifelse(x<=-0.9,1,0)
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
    rownames(cna) <- samples_temp

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
      # filter(rownames(cna) == "P-0000244-T01-IM3") %>%
      mutate_all(~ as.numeric(gsub(" ","",as.character(.)))) %>%
      mutate_all(
        ~ case_when(
          . == 0 ~ "NEUTRAL",
          . %in% c(-1.5,-1) ~ "LOH",
          . == 1 ~ "GAIN",
          . == 2 ~ "AMPLIFICATION",
          . == -2 ~ "DELETION",
        )
      ) %>%
      mutate_all(
        ~factor(.,
                levels = c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION")[
                  which(c("NEUTRAL","DELETION","LOH","GAIN","AMPLIFICATION") %in% .)
                  ])
      )


    colnames(cna) <- paste0(colnames(cna),".cna")
    rownames(cna) <- patients

    cna[is.na(cna)] <- "NEUTRAL"
  }

  return(cna)
}
