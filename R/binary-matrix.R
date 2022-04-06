# Main Binary Matrix Function ------------------

#' Enables creation of a binary matrix from a mutation file with
#' a predefined list of samples (rows are samples and columns are genes)
#' @param samples a character vector specifying which samples should be included in the resulting data frame.
#' Default is NULL is which case all samples with an alteration in the mutation, cna or fusions file will be used. If you specify
#' a vector of samples that contain samples not in any of the passed genomic data frames, 0's (or NAs when appropriate if specifying a panel) will be
#' returned for every column of that patient row.
#' @param mutation A data frame of mutations in the format of a maf file.
#' @param mut_type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param snp_only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include_silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @param fusion A data frame of fusions. If not NULL the outcome will be added to the matrix with columns ending in ".fus".
#' Default is NULL.
#' @param cna A data frame of copy number alterations. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
#' Default is NULL.
#' @param cna_binary A boolean argument specifying if the cna events should be enforced as binary. In which case separate columns for
#' amplifications and deletions will be created.
#' @param cna_relax By default this argument is set to FALSE, where only deep deletions (-2) and amplifications (2) will be annotated as events. When set to TRUE all deletions (-1 shallow and -2 deep) are counted as an event same for all gains (1 gain, 2 amplification) as an event.
#' @param specify_panel a character vector of length 1 with panel id (see gnomeR::gene_panels for available panels), "impact", or "no", or a
#' data frame of `sample_id`-`panel_id` pairs specifying panels for which to insert NAs indicating that gene was not tested.
#' If a single panel id is passed, all genes that are not in that panel (see gnomeR::gene_panels) will be set to NA in results.
#' If `"impact"` passed, impact panel version will be inferred based on each sample_id and NA's annotated accordingly for each sample
#' If specific panel IDs passed, genes not tested in that panel will be should be considered. Default When TRUE NAs will fill the cells for genes
#' Default is "no" which returns data as is with no NA annotation.
#' If you wish to specify different panels for each sample, pass a data frame (with all samples included) with columns: `sample_id`, and `panel_id`. Each sample will be
#' annotated with NAs according to that specific panel.
#' @param rm_empty boolean specifying if alteration columns with no events founds should be removed. Default is FALSE
#' @param recode_aliases boolean specifying if automated gene name alias matching should be done. Default is TRUE. When TRUE
#' the function will check for genes that may have more than 1 name in your data using the aliases im gnomeR::impact_gene_info alias column
#' @param ... Further arguments passed to the oncokb() function such a token
#' @return a data frame with sample_id and alteration binary columns with values of 0/1
#' @export
#' @examples
#' mut.only <- binary_matrix(mutation = mut)
#'
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#'
#' bin.mut <- binary_matrix(samples = samples, mutation = mut,
#' mut_type = "SOMATIC" ,snp_only = FALSE,
#' include_silent = FALSE)
#' bin.mut <- binary_matrix(samples = samples, mutation = mut,
#' mut_type = "SOMATIC", snp_only = FALSE,
#' include_silent = FALSE,
#' cna_relax = TRUE, specify_panel = "no", rm_empty = FALSE)
#' @import dplyr
#' @import dtplyr
#' @import stringr


binary_matrix <- function(samples=NULL,

                          mutation = NULL,
                          mut_type = "SOMATIC",
                          snp_only = FALSE,
                          include_silent = FALSE,

                          fusion = NULL,

                          cna = NULL,
                          cna_binary = TRUE,
                          cna_relax = FALSE,

                          specify_panel = "no",
                          rm_empty = FALSE,
                          recode_aliases = TRUE,
                          ...){

  genie_gene_info <- gnomeR::genie_gene_info
  impact_gene_info <- gnomeR::impact_gene_info
  pathways <- gnomeR::pathways
  gene_panels <- gnomeR::gene_panels

  # Check Arguments ------------------------------------------------------------

  if(is.null(mutation) && is.null(fusion) && is.null(cna)) {
    cli::cli_abort("You must provided at least one of the three following arguments: {.code mutation}, {.code fusion} or {.code cna}.")
  }

  # Check that mutation, fusion, cna is data.frame
  is_df <- purrr::map(
    list(mutation = mutation, fusion = fusion, cna = cna),
    ~dplyr::case_when(
      !is.null(.x) ~ "data.frame" %in% class(.x))
    ) %>%
    purrr::compact()

  not_df <- names(is_df[which(is_df == FALSE)])

  if(length(not_df) > 0) {
    cli::cli_abort("{.code {not_df}} must be a data.frame")
  }


  # if samples not passed we will infer it from data frames
  switch(is.null(samples),
         cli::cli_alert_info("{.code samples} argument is {.code NULL}. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in {.field mutation}, {.field fusion} or {.field cna} data frames"))

  # * Specify_panel must be a known character or data frame with specified column-----

  # make tibbles into data.frames - idk if this is needed, could change switch to ifelse I think a alternative
  if("tbl" %in% class(specify_panel)) {
    specify_panel <- as.data.frame(specify_panel)
  }

  specify_panel <-
    switch(class(specify_panel),
         "character" = {

           choices_arg = c("no", "impact", "IMPACT", gene_panels$gene_panel)
           match.arg(specify_panel, choices = choices_arg)
         },

         "data.frame" = {

           specify_panel %>%
             purrr::when(
               length(setdiff(c(specify_panel$gene_panel), gene_panels$gene_panel)) > 0 ~
                           cli::cli_abort("Panels not known: {.val {setdiff(c(specify_panel$gene_panel), gene_panels$gene_panel)}}. See {.code  gnomeR::gene_panels} for known panels, or skip annotation with {.code specify_panel = 'no'} or indicating {.code 'no'} for those samples in {.field panel_id} column of sample_id-panel_id pair data frame"),
                         TRUE ~ .)},

           cli::cli_abort("{.code specify_panel} must be a character vector of length 1 or a data frame.")
         )


  # * Mutation mutation checks  --------
  mutation <- switch(!is.null(mutation),
                     check_mutation_input(mutation = mutation))

  # * Fusion checks  ----------
  fusion <- switch(!is.null(fusion),
                   check_fusion_input(fusion))

  # * CNA checks  ------------
  cna <- switch(!is.null(cna), {
    check_cna_input(cna)

    # ** This could go in cna check function
    # remove samples with zero events-  this makes CNA match the way mutation/fusion files include samples (which is samples with events only)

    zero_cna <- cna %>%
      purrr::map_at(., vars(-c("Hugo_Symbol")), ~all(.x == 0))

    zero_cna$Hugo_Symbol <- FALSE

    zero_cna <- unlist(zero_cna, use.names = TRUE)
    zero_cna <- zero_cna[zero_cna]
    cna <- cna %>%
      select(-c(names(zero_cna)))

  })

  # ** This assumes regular cna format
  # ** This may not work with research samples!
  cna_samples <- cna %>% purrr::when(!is.null(.) ~ {
                         names(.)[-which(names(.) =='Hugo_Symbol')] %>%
                           str_replace_all(fixed("."), "-")
                         },
                         TRUE ~ NULL)



  #  Make Final Sample List ----------------------------------------------------


  # If user doesn't pass a vector, use samples in files as final sample list
    samples_final <- samples %||%
      c(mutation$Tumor_Sample_Barcode,
        fusion$Tumor_Sample_Barcode,
        cna_samples) %>%
      as.character() %>%
      unique()

    # Binary matrix for each data type ----------------------------------------------
    mutation_binary_df <- switch(!is.null(mutation),
                              .mutations_binary_matrix(mutation = mutation,
                                                       samples = samples_final,
                                                       mut_type = mut_type,
                                                       snp_only = snp_only,
                                                       include_silent = include_silent,
                                                       specify_panel = specify_panel,
                                                       recode_aliases = recode_aliases))


  # fusions
  fusion_binary_df <-  switch(!is.null(fusion),
                           .fusions_binary_matrix(fusion = fusion,
                                                  samples = samples_final,
                                                  specify_panel = specify_panel,
                                                  recode_aliases = recode_aliases))


  # cna
  cna_binary_df <- switch(!is.null(cna),
                       .cna_binary_matrix(cna = cna,
                                          samples = samples_final,
                                          cna_binary = cna_binary,
                                          cna_relax = cna_relax,
                                          specify_panel = specify_panel,
                                          recode_aliases = recode_aliases))

  # put them all together
  all_binary <- bind_cols(list(mutation_binary_df,
                               fusion_binary_df,
                               cna_binary_df))

  # Platform-specific NA Annotation ------

  # we've already checked the arg is valid
  # If character, make into data frame sample-panel pair to input in function
  if(is.character(specify_panel)) {

    sample_panel_pair <- switch(specify_panel,
      "impact" = specify_impact_panels(all_binary),
      "no" = {
        rownames(all_binary) %>%
          as.data.frame() %>%
          stats::setNames("sample_id") %>%
          mutate(panel_id = "no")
      },

      rownames(all_binary) %>%
        as.data.frame() %>%
        stats::setNames("sample_id") %>%
        mutate(panel_id = specify_panel)
    )
    # create data frame of sample IDs

  }

  all_binary <- annotate_any_panel(sample_panel_pair, all_binary)

  # Remove Empty Columns ------
  not_all_na <- which(
      apply(all_binary, 2,
                     function(x){
                       length(unique(x[!is.na(x)]))
                       }) > 1)

  if(rm_empty) {
    all_binary <- all_binary[, which(not_all_na > 1)]
  }

  # Warnings and Attributes --------

  # return omitted zero  samples as warning/attribute

  # samples_omitted <- setdiff(samples, samples_final)
  #
  # if(sum(samples_omitted) > 0) {
  #   attr(all_binary, "omitted") <- samples_omitted
  #
  #   cli::cli_warn("i" = "{length(samples_omitted)} samples were omitted due to  not have any mutations found in the mutation file.
  #                 To view these samples see {.code attr(<your_df>, 'omitted')}")
  # }

  return(all_binary)
}


# Mutations Binary Matrix -----------------------------------------------------

#' Make Binary Matrix From Mutation data frame
#'
#' @inheritParams binary_matrix
#'
#' @return a data frame
#' @export
#'
.mutations_binary_matrix <- function(mutation,
                                     samples,
                                     mut_type,
                                     snp_only,
                                     include_silent,
                                     specify_panel,
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
  mut <- as.data.frame(matrix(0L, nrow=length(samples),
                              ncol=length(unique(mutation$Hugo_Symbol))))

  colnames(mut) <- unique(mutation$Hugo_Symbol)
  rownames(mut) <- samples

  # populate matrix
  for(i in samples){
    genes <- mutation$Hugo_Symbol[mutation$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0) {
      mut[match(i, rownames(mut)), match(unique(as.character(genes)), colnames(mut))] <- 1
      }
  }

  return(mut)
}



# Fusions Binary Matrix -----------------------------------------------------

#' Make Binary Matrix From Fusion data frame
#'
#' @inheritParams binary_matrix
#'
#' @return a data frame
#'
.fusions_binary_matrix <- function(fusion,
                                 samples,
                                 specify_panel,
                                 recode_aliases){


  fusion <- fusion %>%
    filter(.data$Tumor_Sample_Barcode %in% samples)

  # create empty data frame -----
  fusions_out <- as.data.frame(matrix(0L, nrow=length(samples),
                              ncol=length(unique(fusion$Hugo_Symbol))))


  colnames(fusions_out) <- unique(fusion$Hugo_Symbol)
  rownames(fusions_out) <- samples

  # populate matrix
  for(i in samples){
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

#' Make Binary Matrix From CNA data frame
#'
#' @inheritParams binary_matrix
#'
#' @return a data frame
#'
.cna_binary_matrix <- function(cna,
                               samples,
                               cna_binary,
                               cna_relax,
                               specify_panel,
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

  # filter those in final samples list
  cna <- cna[rownames(cna) %in% samples,]

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
    if(length(which(is.na(match(samples,rownames(cna))))) > 0){
      missing <- samples[which(is.na(match(samples,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
    }
    cna <- cna[match(samples,rownames(cna)),]
    cna[is.na(cna)] <- 0
  }

  if(!cna_binary){

    # add missing
    if(length(which(is.na(match(samples,rownames(cna))))) > 0){
      missing <- samples[which(is.na(match(samples,rownames(cna))))]
      add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
      rownames(add )  <- missing
      colnames(add) <- colnames(cna)
      cna <- as.data.frame(rbind(cna,add))
      cna <- cna[match(samples,rownames(cna)),]
    }
    cna <- cna[match(samples,rownames(cna)),]

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
    rownames(cna) <- samples

    cna[is.na(cna)] <- "NEUTRAL"
  }

  return(cna)
}

#binary_matrix(mutation = gnomeR::mut, specify_panel = "no")
