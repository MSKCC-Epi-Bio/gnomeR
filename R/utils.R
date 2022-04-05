#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param mutation Raw maf dataframe containing alteration data
#' @param ... other arguments passed from binary_matrix() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#' check_mutation_input(mutation = gnomeR::mut)
#'
check_mutation_input <- function(mutation, ...)  {

  impact_gene_info <- gnomeR::impact_gene_info
  arguments <- list(...)


  # Check for Fusions-  Old API used to return fusions --------------
  fusions_in_maf <- mutation %>%
    filter(.data$Variant_Classification %in% c("Fusion", "fusion"))

  if(nrow(fusions_in_maf) > 0) {
    cli::cli_abort("It looks like you have fusions in your mutation data frame. These need to be passed to the `fusions` argument. ")
  }

  # Check required columns & data types ------------------------------------------
  required_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  column_names <- colnames(mutation)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Make sure they are character
  mutation <- mutation %>%
    mutate(Tumor_Sample_Barcode = as.character(.data$Tumor_Sample_Barcode),
           Hugo_Symbol = as.character(.data$Hugo_Symbol))

  # * Check suggested columns --------

  # Mutation_Status ---
  if(!("Mutation_Status" %in% column_names)) {
    cli::cli_warn("The following columns are missing in your mutations data: {.field Mutation_Status}. It will be assumed that
            all variants are {.val SOMATIC}.")

    mutation <- mutation %>%
      mutate(Mutation_Status = "SOMATIC")
  }

  # Variant_Type ---
  if(!("Variant_Type" %in% column_names) ) {

    mutation <- mutation %>%
      purrr::when(
        ("Reference_Allele" %in%  column_names) & ("Tumor_Seq_Allele2" %in% column_names) ~

          mutation %>%
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
                              To proceed, add a column specifying {.field Variant_Type} (e.g. {.code mutate(<your-mutation-df>, Variant_Type = 'SNP')}")
      )


    cli::cli_warn("Column {.field Variant_Type} is missing from your data. We inferred variant types using {.field Reference_Allele} and {.field Tumor_Seq_Allele2} columns")

  }
  return(mutation)


}


#' Check fusion data frame to ensure columns are correct
#'
#' @param fusion a fusion data frame
#' @param ... other arguments passed from binary_matrix()
#'
#' @return a checked data frame
#' @export
#'
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
    mutate(Tumor_Sample_Barcode = as.character(.data$Tumor_Sample_Barcode),
           Hugo_Symbol = as.character(.data$Hugo_Symbol))

  return(fusion)
}



#' Check CNA data frame to ensure columns are correct
#'
#' @param cna a cna data frame
#' @param ... other arguments passed from binary_matrix()
#'
#' @return a checked data frame
#' @export
#'
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
    mutate(Hugo_Symbol = as.character(.data$Hugo_Symbol))


  return(cna)
}

# Recode Gene Aliases---------------------

#' Recode Hugo Symbol Column
#'
#' @param genomic_df a binary_matrix object
#' @param ... Other things passed
#'
#' @return A dataframe with a recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' recode_alias(gnomeR::mut)
#'
recode_alias <- function(genomic_df, ...) {

  # get table of gene aliases (internal data)
    alias_table <- tidyr::unnest(gnomeR::impact_gene_info, cols = .data$alias) %>%
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
resolve_alias <- function(gene_to_check, alias_table) {

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


#' Reformat CNA from maf version to wide version
#'
#' @param cna a cna dataframe in maf format
#' @param patients a list of patients to include
#'
#' @return a dataframe of reformatted CNA alteration
#' @export
#'
reformat_cna <- function(cna, patients = patients) {

  # recreate original format
  temp <- as.data.frame(matrix(0L ,ncol = length(patients)+1,
                               nrow = length(unique(cna$Hugo_Symbol))))

  colnames(temp) <- c("Hugo_Symbol",patients)
  temp[,1] <- unique(cna$Hugo_Symbol)

  for(i in patients){
    temp[match(as.character(unlist(cna %>%
                                     filter(.data$sampleId %in% i) %>%
                                     select(.data$Hugo_Symbol))),temp[,1]),
         match(i, colnames(temp))] <- as.numeric(unlist(cna %>%
                                                          filter(.data$sampleId %in% i) %>%
                                                          select(.data$alteration)))
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



#' IMPACT Panel Annotation of NA's
#'
#' @param binary_matrix a processed binary_matrix
#'
#' @return  a data frame iwth NAs inserted for genes not tested for given panel versions
#' @export
#'
#'
annotate_impact_missing <- function(binary_matrix) {

  gene_panels <- gnomeR::gene_panels

  # create data frame of sample IDs
  sample_panel_pair <- rownames(binary_matrix) %>%
    as.data.frame() %>%
    stats::setNames("sample_id")

  any_impact <- sum(stringr::str_detect(sample_panel_pair$sample_id,
                                        "-IM|-IH"))

  # check if there are any impact (based on sample ID)
  switch(any_impact == 0,
         cli::cli_abort("There are no IMPACT samples recognized (based on sample_id). If you wish to annotate
                   missingness based on a panel, please pass a data.frame of sample_ids and corresponding panels."))


  # get which IMPACT panel
  sample_panel_pair <- sample_panel_pair %>%
    mutate(gene_panel = case_when(
      stringr::str_detect(.data$sample_id, "-IM3") ~ "IMPACT341",
      stringr::str_detect(.data$sample_id, "-IM5") ~ "IMPACT410",
      stringr::str_detect(.data$sample_id, "-IM6") ~ "IMPACT468",
      stringr::str_detect(.data$sample_id, "-IM7") ~ "IMPACT505",
      TRUE ~ "none"
    ))

  sample_panel_pair_nest <- sample_panel_pair %>%
    group_by(.data$gene_panel) %>%
    summarise(samples_in_panel = list(.data$sample_id))

  # pull genes for given panels
  panels_needed <- unique(sample_panel_pair_nest$gene_panel)

  # has sample IDs and genes for each panel
  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    left_join(gene_panels) %>%
    select(-.data$entrez_ids_in_panel)

  user_data_genes <- gsub(".fus|.Del|.Amp|.cna", "", colnames(binary_matrix))

  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    mutate(na_genes_raw = purrr::map(.data$genes_in_panel,
                                 ~unique(setdiff(user_data_genes, .x)))) %>%
    mutate(na_genes = purrr::map(.data$na_genes_raw,
                                     ~c(
                                       .x,
                                       paste0(.x, ".fus"),
                                       paste0(.x, ".Del"),
                                       paste0(.x, ".Amp"),
                                       paste0(.x, ".cna")
                                     )))


  annotated_data <- purrr::pmap_df(sample_panel_pair_nest,
                                   annotate_panel,
                                   binary_matrix = binary_matrix)

  return(annotated_data)
}


#' Utility function  to insert NA's According to Panel
#'
#' @param binary_matrix a processed binary matrix
#' @param gene_panel name of gene panel
#' @param samples_in_panel samples to be annotated for each panel
#' @param na_genes genes to make NA
#' @param ... other args passed
#'
#' @return an annotated data frame
#' @export
#'
#'
annotate_panel <- function(binary_matrix,
                           gene_panel,
                           samples_in_panel,
                           na_genes, ...) {

  mut_sub <- binary_matrix[samples_in_panel, ]
  mut_sub[,stats::na.omit(match(na_genes, colnames(mut_sub)))] <- NA

  return(mut_sub)

}





