# Main Binary Matrix Function ------------------

#' binmat
#' \%lifecycle{experimental}
#' Enables creation of a binary matrix from a mutation file with
#' a predefined list of patients (rows are patients and columns are genes)
#' @param patients a character vector specifying which patients should be included in the resulting data frame.
#' Default is NULL is which case all patients within the mutation, cna or fusions file will be used. If you specify
#' a vector of patients that contain patients not in any of the passed genomic data frames, 0's (or NAs when appropriate if specifying a panel) will be
#' returned for every column of that patient row.
#' @param mutation A data frame of mutations in the format of a maf file.
#' @param mut_type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param snp_only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include_silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @param fusion An optional data frame of fusions. If not NULL the outcome will be added to the matrix with columns ending in ".fus".
#' Default is NULL.
#' @param cna An optional CNA data frame. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
#' Default is NULL.
#' @param cna_binary A boolean argument specifying if the cna events should be enforced as binary. In which case separate columns for
#' amplifications and deletions will be created.
#' @param cna_relax By default this argument is set to FALSE, where only deep deletions (-2) and amplifications (2) will be annotated as events.
#'  When set to FTRUE all deletions (-1 shallow and -2 deep) are counted as an event same for all gains (1 gain, 2 amplification) as an event.
#' @param specify_panel boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is TRUE.
#' @param rm_empty boolean specifying if columns with no events founds should be removed. Default is TRUE.
#' @param recode_aliases bolean specifying if automated gene name alias matching should be done. Default is TRUE. When TRUE
#' the function will check for genes that may have more than 1 name in your data using the aliases im gnomeR::impact_gene_info alias column
#' @param col.names character vector of the necessary columns to be used. By default: col.names = c(Tumor_Sample_Barcode = NULL,
#'  Hugo_Symbol = NULL, Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL)
#' @param oncokb boolean specfiying if mutation file should be oncokb annotated. Default is FALSE.
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
#' @examples
#' # mut.only <- binmat(mutation = mut)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,mutation = mut,
#' mut_type = "SOMATIC",snp_only = FALSE,
#' include_silent = FALSE, specify_platform = FALSE)
#' bin.mut <- binmat(patients = patients, mutation = mut,
#' mut_type = "SOMATIC",snp_only = FALSE,
#' include_silent = FALSE,
#' cna_relax = TRUE, specify_platform = FALSE,
#'  set.plat = "410", rm_empty = FALSE)
#' @import dplyr
#' @import dtplyr
#' @import stringr


binary_matrix <- function(patients=NULL,

                          mutation = NULL,
                          mut_type = "SOMATIC",

                          snp_only = FALSE,
                          include_silent = FALSE,

                          fusion = NULL,
                          cna = NULL,
                          cna_binary = TRUE,
                          cna_relax = FALSE,

                          specify_panel = "impact",
                          rm_empty = TRUE,
                          pathway = FALSE,
                          recode_aliases = TRUE,

                   col.names = c(Tumor_Sample_Barcode = NULL, Hugo_Symbol = NULL,
                                 Variant_Classification = NULL, Mutation_Status = NULL,
                                 Variant_Type = NULL),

                   oncokb = FALSE,
                   keep_onco = c("Oncogenic","Likely Oncogenic","Predicted Oncogenic"),
                   token = '',
                   sample_panels = NULL,...){

  genie_gene_info <- gnomeR::genie_gene_info
  impact_gene_info <- gnomeR::impact_gene_info
  panel_names <- gnomeR::panel_names
  pathways <- gnomeR::pathways

  # Check Arguments ------------------------------------------------------------

  if(is.null(mutation) && is.null(fusion) && is.null(cna)) {
    cli::cli_abort("You must provided at least one of the three following arguments: {.code mutation}, {.code fusion} or {.code cna}.")
  }

  # Check that mutation, fusion, cna is data.frame
  is_df <- purrr::map(
    list(mutation = mutation, fusion = fusion, cna= cna),
    ~dplyr::case_when(
      !is.null(.x) ~ "data.frame" %in% class(.x))
    ) %>%
    purrr::compact()

  not_df <- names(is_df[which(is_df == FALSE)])

  if(length(not_df) > 0) {
    cli::cli_abort("{.code {not_df}} must be a data.frame")
  }


  # * Mutation mutation checks  --------
  mutation <- switch(!is.null(mutation),
                     check_mutation_input(maf = mutation))

  # * Fusion checks  ----------
  fusion <- switch(!is.null(fusion),
                   check_fusion_input(fusion))

  # * CNA checks  ------------
  cna <- switch(!is.null(cna),
                check_cna_input(cna))


  #  Make Final Sample List ----------------------------------------------------

  # If user doesn't pass a vector, use samples in files as final sample list
    patients_final <- patients %||%
      c(mutation$Tumor_Sample_Barcode,
        fusion$Tumor_Sample_Barcode,
        cna$Tumor_Sample_Barcode) %>%
      as.character() %>%
      unique()

    # Binary matrix for each data type ----------------------------------------------
    mutation_binary_df <- switch(!is.null(mutation),
                              .mutations_binary_matrix(mutation = mutation,
                                                       patients = patients_final,
                                                       mut_type = mut_type,
                                                       snp_only = snp_only,
                                                       include_silent = include_silent,
                                                       specify_platform = specify_platform,
                                                       recode_aliases = recode_aliases))

  # fusions
  fusion_binary_df <-  switch(!is.null(fusion),
                           .fusions_binary_matrix(fusion = fusion,
                                                  patients = patients_final,
                                                  specify_platform = specify_platform,
                                                  recode_aliases = recode_aliases))


  # cna
  cna_binary_df <- switch(!is.null(cna),
                       .cna_binary_matrix(cna = cna,
                                          patients = patients_final,
                                          cna_binary = cna_binary,
                                          cna_relax = cna_relax,
                                          specify_platform = specify_platform,
                                          recode_aliases = recode_aliases))


  all_binary <- bind_cols(list(mutation_binary_df,
                               fusion_binary_df,
                               cna_binary_df))

  # Platform-specific Missingness Annotation ------

  all_binary <- switch(specify_panel,
                       "impact" = annotate_impact_missing(all_binary),
                       "no" = all_binary)



  return(all_binary)

  # Remove Empty Columns ------
  not_all_na <- which(
      apply(all_binary, 2,
                     function(x){
                       length(unique(x[!is.na(x)]))
                       })>1)

  if(rm_empty)

    mut <- mut[,which(apply(all_binary,2,function(x){length(unique(x[!is.na(x)]))})>1)]
}


# Mutations Binary Matrix -----------------------------------------------------

.mutations_binary_matrix <- function(mutation,
                                     patients = patients_final,
                                     mut_type,
                                     snp_only,
                                     include_silent,
                                     specify_platform,
                                     recode_aliases = recode_aliases){

  # apply filters --------------
  mutation <- mutation %>%
    purrr::when(
      snp_only ~ filter(., .data$Variant_Type == "SNP"),
      ~ .) %>%

    purrr::when(
      !include_silent ~ filter(., .data$Variant_Classification != "Silent"),
      ~ .) %>%

    purrr::when(
      tolower(mut_type) == "all" ~ .,
      TRUE ~ filter(., .data$Mutation_Status %in% mut_type),
      ~ .)


  # create empty data.frame to hold results -----
  mut <- as.data.frame(matrix(0L, nrow=length(patients),
                              ncol=length(unique(mutation$Hugo_Symbol))))

  colnames(mut) <- unique(mutation$Hugo_Symbol)
  rownames(mut) <- patients

  # populate matrix
  for(i in patients){
    genes <- mutation$Hugo_Symbol[mutation$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0) {
      mut[match(i, rownames(mut)), match(unique(as.character(genes)), colnames(mut))] <- 1
      }
  }


  zero_mutations <- apply(mut, 1, function(x){ sum(x)== 0 })
  zero_mutations <- zero_mutations[zero_mutations == TRUE]


  # return omitted zero mutations samples as warning/attribute
  if(sum(zero_mutations) > 0) {
    attr(mut, "omitted") <- zero_mutations

    cli::cli_warn("{length(zero_mutations)} samples did not have any mutations found in the mutation file. To view these samples see {.code attr(mut, 'omitted')}")
 }
  return(mut)
}



# Fusions Binary Matrix -----------------------------------------------------

.fusions_binary_matrix <- function(fusion,
                                 patients = patients_final,
                                 specify_platform,
                                 recode_aliases){


  fusion <- as_tibble(fusion) %>%
    filter(.data$Tumor_Sample_Barcode %in% patients)

  # create empty data frame -----
  fusions_out <- as.data.frame(matrix(0L, nrow=length(patients),
                              ncol=length(unique(fusion$Hugo_Symbol))))


  colnames(fusions_out) <- unique(fusion$Hugo_Symbol)
  rownames(fusions_out) <- patients

  # populate matrix
  for(i in patients){
    genes <- fusion$Hugo_Symbol[fusion$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0)
      fusions_out[match(i,rownames(fusions_out)),
                 match(unique(as.character(genes)),colnames(fusions_out))] <- 1

  }

  # add .fus suffix on columns
  colnames(fusions_out) <- paste0(colnames(fusions_out),".fus")
  return(fusions_out)
}


# CNA Binary Matrix -----------------------------------------------------

.cna_binary_matrix <- function(cna,
                               patients = patients_final,
                               cna_binary,
                               cna_relax,
                               specify_platform,
                               recode_aliases){



  # If more than 1 row per gene, combine rows
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

  # flip
  cna <- as.data.frame(t(cna))

  # fix names
  rownames(cna) <- gsub("\\.","-",rownames(cna))

  # filter those in final patients list
  cna <- cna[rownames(cna) %in% patients,]

  # If cna binary
  samples_temp <- rownames(cna)

  if(cna_binary){
    temp <- do.call("cbind",apply(cna,2,function(x){
      if(cna_relax){
        yA <- ifelse(as.numeric(x)>=0.9,1,0)
        yD <- ifelse(as.numeric(x)<=-0.9,1,0)
      }
      if(!cna_relax){
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

  if(!cna_binary){

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


    # add cna annotation
    colnames(cna) <- paste0(colnames(cna),".cna")
    rownames(cna) <- patients

    cna[is.na(cna)] <- "NEUTRAL"
  }

  return(cna)
}

#binary_matrix(mutation = gnomeR::mut, specify_panel = "no")
