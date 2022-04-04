#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param ... Further arguments parsed through binmat() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#' check_mutation_input(mut,recode.aliases = TRUE)
#'
check_mutation_input <- function(maf, ...)  {

  impact_gene_info <- gnomeR::impact_gene_info
  arguments <- list(...)


  # Check for Fusions-  Old API used to return fusions --------------
  fusions_in_maf <- maf %>%
    filter(.data$Variant_Classification %in% c("Fusion", "fusion"))

  if(nrow(fusions_in_maf) > 0) {
    cli::cli_abort("It looks like you have fusions in your maf. These need to be passed to the `fusions` argument. ")
  }

  # Check required columns & data types ------------------------------------------
  required_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  column_names <- colnames(maf)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Make sure they are character
  maf <- maf %>%
    mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
           Hugo_Symbol = as.character(Hugo_Symbol))

  # * Check suggested columns --------

  # Mutation_Status ---
  if(!("Mutation_Status" %in% column_names)) {
    cli::cli_warn("The following columns are missing in your mutations data: {.field Mutation_Status}. It will be assumed that
            all variants are {.val SOMATIC}.")

    maf <- maf %>%
      mutate(Mutation_Status = "SOMATIC")
  }

  # Variant_Type ---
  if(!("Variant_Type" %in% column_names) ) {

    maf <- maf %>%
      purrr::when(
        ("Reference_Allele" %in%  column_names) & ("Tumor_Seq_Allele2" %in% column_names) ~

          maf %>%
            mutate(
              Reference_Allele = as.character(.data$Reference_Allele),
              Tumor_Seq_Allele2 = as.character(.data$Tumor_Seq_Allele2),
              Variant_Type = case_when(
                .data$Reference_Allele %in% c("A","T","C","G") &
                .data$Tumor_Seq_Allele2 %in% c("A","T","C","G") ~ "SNP",
                nchar(.data$Tumor_Seq_Allele2) < nchar(.data$Reference_Allele) |
                .data$Tumor_Seq_Allele2 == "-" ~ "DEL",
                .data$Reference_Allele == "-" |
                nchar(.data$Tumor_Seq_Allele2) > nchar(.data$Reference_Allele) ~ "INS",
                nchar(.data$Reference_Allele) == 2 & nchar(.data$Tumor_Seq_Allele2) == 2 ~ "DNP",
                nchar(.data$Reference_Allele) == 3 & nchar(.data$Tumor_Seq_Allele2) == 3 ~ "TNP",
                nchar(.data$Reference_Allele) > 3 & nchar(.data$Tumor_Seq_Allele2) == nchar(.data$Reference_Allele) ~ "ONP",
                TRUE ~ "Undefined")),

        TRUE ~ cli::cli_abort("Column {.field Variant_Type} is missing from your data and {.field Reference_Allele} and {.field Tumor_Seq_Allele2}
                              columns were not available from which to infer variant type.
                              To proceed, add a column specifying {.field Variant_Type} (e.g. {.code mutate(<your-maf>, Variant_Type = 'SNP')}")
      )


    cli::cli_warn("Column {.field Variant_Type} is missing from your data. We inferred variant types using {.field Reference_Allele} and {.field Tumor_Seq_Allele2} columns")

  }
  return(maf)


}


check_fusion_input <- function(fusion, ...)  {

  impact_gene_info <- gnomeR::impact_gene_info
  arguments <- list(...)


  # Check required columns & data types ------------------------------------------
  required_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol")
  column_names <- colnames(fusion)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Make sure they are character
  fusion <- fusion %>%
    mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
           Hugo_Symbol = as.character(Hugo_Symbol))


}


#' Check CNA Data
#'
#' @param cna
#' @param ...
#'
#' @return a data frame
#' @export
#'
#' @examples
check_cna_input <- function(cna, ...)  {

  impact_gene_info <- gnomeR::impact_gene_info
  arguments <- list(...)

  # Check required columns & data types ------------------------------------------
  required_cols <- c("Hugo_Symbol")
  column_names <- colnames(cna)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Check format of CNA
  api_cols <- c("sampleId", "studyId", "patientId", "alteration")
  in_data <- api_cols[api_cols %in% column_names]

  if(length(in_data) > 0) {
    cli::cli_abort("The following columns are not allowed in your cna data.frame {.field {in_data}}.
                   Do you need to reformat? See {.code ?reformat_cna()}")
  }

  # Make sure Hugo is character
  cna <- cna %>%
    mutate(Hugo_Symbol = as.character(Hugo_Symbol))


}

# Recode Gene Aliases---------------------

#' Recode Hugo Symbol Column
#'
#' @param genomic_df
#'
#' @return A dataframe with a recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' recode_alias(gnomeR::mut)
#'
recode_alias <- function(genomic_df, ...) {

  # get table of gene aliases (internal data)
    alias_table <- tidyr::unnest(impact_gene_info, cols = .data$alias) %>%
      dplyr::select(.data$hugo_symbol, .data$alias)

    # recode aliases
    genomic_df$Hugo_Symbol_Old <- genomic_df$Hugo_Symbol
    genomic_df$Hugo_Symbol <- purrr::map_chr(genomic_df$Hugo_Symbol, ~resolve_alias(gene_to_check = .x,
                                                                      alias_table = alias_table))

    message <- genomic_df %>%
      dplyr::filter(.data$Hugo_Symbol_Old != .data$Hugo_Symbol) %>%
      dplyr::select(.data$Hugo_Symbol_Old, .data$Hugo_Symbol) %>%
      dplyr::distinct()


    if(nrow(message) > 0) {
      vec_recode <- purrr::map2_chr(message$Hugo_Symbol_Old,
                                 message$Hugo_Symbol,
                                 ~paste0(.x, " recoded to ", .y))

      names(vec_recode) <- rep("!", times = length(vec_recode))

      cli::cli_warn(c(
      "To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your dataframe have been recoded (you can supress this with {.code recode.aliases = FALSE}):",
      vec_recode))

  }

  return(genomic_df)
}




#' Utility Function to Extract SNV
#'
#' @param x string
#' @param n number of characters from right
#'
#' @return string
#' @export
#'
#' @examples
#' substrRight("Hello", 2)
#'
substrRight <- function(x, n) {
  x <- as.character(x)
  substr(x, nchar(x) - n + 1, nchar(x))
}



#' Resolve Hugo Symbol Names with Aliases
#'
#' @param gene_to_check hugo_symbol to be check
#' @param alias_table table containing all the aliases
#'
#' @return if the accepted hugo symbol is input, it is returned back.
#' If an alias name is provided, the more common name/more up to date name is returned
#' @export
#'
#' @examples
#' resolve_alias("KMT2D", alias_table = tidyr::unnest(impact_gene_info, cols = alias))
#'
resolve_alias <- function(gene_to_check, alias_table = all_alias_table) {

  if(gene_to_check %in% alias_table$alias) {

    alias_table %>%
      filter(.data$alias == gene_to_check) %>%
      pull(.data$hugo_symbol) %>%
      first() %>%
      as.character()

  } else {
    as.character(gene_to_check)
  }
}


reformat_cna <- function(cna, patients = patients_final) {

  # recreate orginal format
  temp <- as.data.frame(matrix(0L ,ncol = length(patients)+1,
                               nrow = length(unique(cna$Hugo_Symbol))))

  colnames(temp) <- c("Hugo_Symbol",patients)
  temp[,1] <- unique(cna$Hugo_Symbol)

  for(i in patients){
    temp[match(as.character(unlist(cna %>%
                                     filter(.data$sampleId %in% i) %>%
                                     select(.data$Hugo_Symbol))),temp[,1]),
         match(i, colnames(temp))] <- as.numeric(unlist(cna %>% filter(.data$sampleId %in% i) %>% select(.data$alteration)))
  }

  cna <- temp
  rownames(cna) <- cna[,1]
  cna <- cna[,-1]
  cna <- as.data.frame(t(cna))
  cna <- cna[rownames(cna) %in% patients,]




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


}

check_cna_input <- function(cna) {

  # Is it API maf style or traditional with samples as colnames?


}




annotate_impact_missing <- function(binary_matrix) {

  any_impact <- sum(stringr::str_detect(sample_panel_pair$sample_id,
                                        "-IM|-IH"))

  # check if there are any impact (based on sample ID)
  switch(any_impact == 0,
         cli::cli_abort("There are no IMPACT samples recognized (based on sample_id). If you wish to annotate
                   missingness based on a panel, please pass a data.frame of sample_ids and corresponding panels."))

  # create data frame of sample IDs
  sample_panel_pair <- rownames(binary_matrix) %>%
    as.data.frame() %>%
    stats::setNames("sample_id")


  # get which IMPACT panel
  sample_panel_pair <- sample_panel_pair %>%
    mutate(gene_panel = case_when(
      stringr::str_detect(sample_id, "-IM3") ~ "IMPACT341",
      stringr::str_detect(sample_id, "-IM5") ~ "IMPACT410",
      stringr::str_detect(sample_id, "-IM6") ~ "IMPACT468",
      stringr::str_detect(sample_id, "-IM7") ~ "IMPACT505",
      TRUE ~ "none"
    ))

  panel_metadata <- sample_panel_pair %>%
    group_by(gene_panel) %>%
    summarise(samples_in_panel = list(sample_id))

  # pull genes for given panels
  panels <- unique(sample_panel_pair$gene_panel)
  panels <- get_gene_panel(panels) %>%
    transmute(gene_panel = genePanelId,
              hugo_symbol = hugoGeneSymbol)

  panels <- panels %>%
    group_by(gene_panel) %>%
    summarise(genes_in_panel = list(hugo_symbol))


  # has sample IDs and genes for each panel
  panel_metadata <- panel_metadata %>%
    left_join(panels)


  user_data_genes <- gsub(".fus|.Del|.Amp|.cna", "", colnames(binary_matrix))

  panel_metadata <- panel_metadata %>%
    mutate(na_genes_raw = purrr::map(genes_in_panel,
                                 ~unique(setdiff(user_data_genes, .x)))) %>%
    mutate(na_genes = purrr::map(na_genes_raw,
                                     ~c(
                                       .x,
                                       paste0(.x, ".fus"),
                                       paste0(.x, ".Del"),
                                       paste0(.x, ".Amp"),
                                       paste0(.x, ".cna")
                                     )))


  annotated_data <- purrr::pmap_df(panel_metadata, annotate_panel,
                                   binary_matrix = all_binary)

  return(annotated_data)
}


annotate_panel <- function(binary_matrix,
                           gene_panel,
                           samples_in_panel,
                           na_genes, ...) {

  mut_sub <- binary_matrix[samples_in_panel, ]
  mut_sub[,stats::na.omit(match(na_genes, colnames(mut_sub)))] <- NA

  return(mut_sub)

}





