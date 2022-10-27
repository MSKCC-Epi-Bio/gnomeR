#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param mutation Raw maf dataframe containing alteration data
#' @param ... other arguments passed from create_gene_binary() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#' sanitize_mutation_input(mutation = gnomeR::mutations)
#'
sanitize_mutation_input <- function(mutation, ...)  {

  arguments <- list(...)

  mutation <- rename_columns(mutation)

  # Check required columns & data types ------------------------------------------
  required_cols <- c("sample_id", "hugo_symbol")
  column_names <- colnames(mutation)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Make sure they are character
  mutation <- mutation %>%
    mutate(sample_id = as.character(.data$sample_id),
           hugo_symbol = as.character(.data$hugo_symbol))

  # Check for Fusions-  Old API used to return fusions --------------
  fusions_in_maf <- mutation %>%
    filter(.data$variant_classification %in% c("Fusion", "fusion"))

  if(nrow(fusions_in_maf) > 0) {
    cli::cli_abort("It looks like you have fusions in your mutation data frame. These need to be passed to the `fusions` argument. ")
  }

  # * Check suggested columns --------

  # Mutation_Status ---
  if(!("mutation_status" %in% column_names)) {
    cli::cli_warn("A {.field mutation_status} column was not found. It will be assumed that
            all variants are {.val SOMATIC}, or check your data follows naming guidelines in {.code gnomer::names_df}")

    mutation <- mutation %>%
      mutate(mutation_status = "SOMATIC")
  }

  # Variant_Type ---
  if(!("variant_type" %in% column_names) ) {

    mutation <- mutation %>%
      purrr::when(
        ("reference_allele" %in%  column_names) & ("tumor_seq_allele2" %in% column_names) ~

          mutation %>%
            mutate(
              reference_allele = as.character(.data$reference_allele),
              tumor_seq_allele2 = as.character(.data$tumor_seq_allele2),
              variant_type = case_when(
                .data$reference_allele %in% c("A","T","C","G") &
                .data$tumor_seq_allele2 %in% c("A","T","C","G") ~ "SNP",
                nchar(.data$tumor_seq_allele2) < nchar(.data$reference_allele) |
                .data$tumor_seq_allele2 == "-" ~ "DEL",
                .data$reference_allele == "-" |
                nchar(.data$tumor_seq_allele2) > nchar(.data$reference_allele) ~ "INS",
                nchar(.data$reference_allele) == 2 & nchar(.data$tumor_seq_allele2) == 2 ~ "DNP",
                nchar(.data$reference_allele) == 3 & nchar(.data$tumor_seq_allele2) == 3 ~ "TNP",
                nchar(.data$reference_allele) > 3 & nchar(.data$tumor_seq_allele2) == nchar(.data$reference_allele) ~ "ONP",
                TRUE ~ "Undefined")),

        TRUE ~ cli::cli_abort("Column {.field variant_type} is missing from your data and {.field reference_allele} and {.field tumor_seq_allele2}
                              columns were not available from which to infer variant type.
                              To proceed, add a column specifying {.field variant_type} (e.g. {.code mutate(<your-mutation-df>, variant_type = 'SNP')}")
      )


    cli::cli_warn("Column {.field variant_type} is missing from your data. We inferred variant types using {.field reference_allele} and {.field tumor_seq_allele2} columns")

  }
  return(mutation)


}


#' Check fusion data frame to ensure columns are correct
#'
#' @param fusion a fusion data frame
#' @param ... other arguments passed from create_gene_binary()
#'
#' @return a checked data frame
#' @export
#' @examples
#' fus <- sanitize_fusion_input(fusion = gnomeR::sv)
#'
sanitize_fusion_input <- function(fusion, ...)  {

  arguments <- list(...)

  fusion <- rename_columns(fusion)

  # Check required columns & data types ------------------------------------------

  column_names <- colnames(fusion)
  # check for hugo symbol

  if(!("sample_id" %in% column_names) > 0) {
    cli::cli_abort("No sample ID column found.")
  }

  if(!("hugo_symbol" %in% column_names |
       "site_1_entrez_gene_id" %in% column_names)) {

    cli::cli_abort("No hugo symbol column found. See `gnomeR::names_df` for accepted column names")
  }


  # Make sure they are character
  fusion <- fusion %>%
    mutate(sample_id = as.character(.data$sample_id)) %>%
    purrr::when(
      "hugo_symbol" %in% column_names ~
        mutate(., "hugo_symbol" = as.character(.data$hugo_symbol)),
      TRUE ~
        mutate(., "site_1_entrez_gene_id" = as.character(.data$site_1_entrez_gene_id)))

  return(fusion)
}


#' Internal function to recode numeric CNA alteration values to factor values
#'
#' @param cna a maf (long) form data set of CNAs. Must include an alteration column.
#'
#' @return a recoded CNA data set with factor alteration values
#'
#'


.recode_cna_alterations <- function(cna){


  if(!("alteration" %in% colnames(cna))) {
    cli::cli_abort("An alteration column is missing from your cna data. Use pivot_cna_longer(),
                   which will recode alterations, instead if dataset is in API format.")
  }

  # Make sure hugo & alteration is character
  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(str_trim(as.character(.data$alteration))))


  #assess levels of alteration
  levels_in_data <- names(table(cna$alteration))

  allowed_chr_levels <- c(
    "neutral" = "0",
    "homozygous deletion" = "-2",
    "loh" = "-1.5", #this is a placeholder until methods cleared up
    "hemizygous deletion" = "-1",
    "gain" = "1",
    "high level amplification" = "2"
  )

  #pull any numbers not in allowed list out
  all_allowed <- c(allowed_chr_levels, names(allowed_chr_levels))
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  #abort if unknown values exist
  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }

  # recode the alteration varaible as factor with those levels
  # and suppress warnings on this
  suppressWarnings(
    cna <- cna %>%
      mutate(alteration = forcats::fct_recode(.data$alteration, !!!allowed_chr_levels))
  )

  return(cna)
}


#' Check CNA data frame to ensure columns are correct
#'
#' @param cna a cna data frame
#' @param ... other arguments passed from create_gene_binary()
#'
#' @return a checked data frame
#' @export
#' @examples
#'
#' cna <- sanitize_cna_input(cna = cna)
#'
sanitize_cna_input <- function(cna, ...)  {

  arguments <- list(...)

  cna <- rename_columns(cna)

  # Check required columns & data types ------------------------------------------
  required_cols <- c("hugo_symbol", "sample_id", "alteration")
  column_names <- colnames(cna)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}.
                   Is your data in wide format? If so, it must be long format. See {.code gnomeR::pivot_cna_long()} to reformat")
  }

  cna <- switch(!is.null(cna), .recode_cna_alterations(cna))


  return(cna)
}

#' Rename columns from API results to work with gnomeR functions
#'
#' @param df_to_check a data frame to check and recode names as needed
#'
#' @return a renamed data frame
#' @export
#' @examples
#'
#' rename_columns(df_to_check = gnomeR::mutations)
#' rename_columns(df_to_check = gnomeR::sv)
#'
rename_columns <- function(df_to_check) {

  names_df_long <- gnomeR::names_df %>%
    select(contains("_column_name")) %>%
    tidyr::pivot_longer(-"internal_column_name")


  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector
  temp_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(temp_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(temp_dict),
                               temp_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check %>%
      dplyr::rename(!!temp_dict)
  }
}




# Recode Gene Aliases---------------------

#' Recode Hugo Symbol Column
#'
#' @param genomic_df a gene_binary object
#' @param ... Other things passed
#'
#' @return A dataframe with a recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' mut <- rename_columns(gnomeR::mutations)
#'
#' colnames(mut)
#'
#' colnames(recode_alias(genomic_df = mut))
#'
recode_alias <- function(genomic_df, ...) {

  # get table of gene aliases (internal data)
    alias_table <- gnomeR::impact_alias_table %>%
      dplyr::select("hugo_symbol", "alias")

    # recode aliases
    genomic_df$hugo_symbol_old <- genomic_df$hugo_symbol
    genomic_df$hugo_symbol <- purrr::map_chr(genomic_df$hugo_symbol,
                                             ~resolve_alias(gene_to_check = .x,
                                                            alias_table = alias_table))

    message <- genomic_df %>%
      dplyr::filter(.data$hugo_symbol_old != .data$hugo_symbol) %>%
      dplyr::select("hugo_symbol_old", "hugo_symbol") %>%
      dplyr::distinct()


    if(nrow(message) > 0) {
      vec_recode <- purrr::map2_chr(message$hugo_symbol_old,
                                 message$hugo_symbol,
                                 ~paste0(.x, " recoded to ", .y))

      names(vec_recode) <- rep("!", times = length(vec_recode))

      cli::cli_warn(c(
      "To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your dataframe have been recoded (you can prevent this with {.code recode_aliases = FALSE}):",
      vec_recode))

  }

    genomic_df <- genomic_df %>%
      select(-"hugo_symbol_old")

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
#' resolve_alias("MLL4", alias_table = impact_alias_table)
#'
resolve_alias <- function(gene_to_check, alias_table) {

  if(gene_to_check %in% alias_table$alias) {

    alias_table %>%
      filter(.data$alias == gene_to_check) %>%
      pull("hugo_symbol") %>%
      first() %>%
      as.character()

  } else {
    as.character(gene_to_check)
  }
}


#' IMPACT Panel Annotation of NA's
#'
#' @param gene_binary a processed gene_binary
#'
#' @return  a data frame iwth NAs inserted for genes not tested for given panel versions
#' @export
#'
#'
specify_impact_panels <- function(gene_binary) {

  gene_panels <- gnomeR::gene_panels

  # create data frame of sample IDs
  sample_panel_pair <- rownames(gene_binary) %>%
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
    mutate(panel_id = case_when(
      stringr::str_detect(.data$sample_id, "-IM3") ~ "IMPACT341",
      stringr::str_detect(.data$sample_id, "-IM5") ~ "IMPACT410",
      stringr::str_detect(.data$sample_id, "-IM6") ~ "IMPACT468",
      stringr::str_detect(.data$sample_id, "-IM7") ~ "IMPACT505",
      TRUE ~ "no"
    ))

  # couldn't detect panel
  unk_impact_panel <- sample_panel_pair %>%
    filter(.data$panel_id == "no")

  if(nrow(unk_impact_panel) > 0) {
    cli::cli_alert("Couldn't infer IMPACT panel version from these sample_ids, therefore no NA panel annotation will be done for these: {unk_impact_panel$sample_id}")
  }

  return(sample_panel_pair)
}

#' Annotate Missing Gene Values According to Specific Panels
#'
#' @param sample_panel_pair a data frame of `sample_id`-`panel_id` pairs specifying panels to use for annotation of each sample
#' @param gene_binary a binary matrix of 0/1 indicating alteration yes/no for each sample
#' @return a gene_binary annotated for missingness
#' @export

annotate_any_panel <- function(sample_panel_pair, gene_binary) {

  # if all "no", leave function
  switch(all(sample_panel_pair$panel_id == "no"),
         return(gene_binary))

  sample_panel_pair_nest <- sample_panel_pair %>%
    group_by(.data$panel_id) %>%
    summarise(samples_in_panel = list(.data$sample_id))

  # pull genes for given panels
  panels_needed <- unique(sample_panel_pair_nest$panel_id)

  # has sample IDs and genes for each panel
  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    left_join(gnomeR::gene_panels, by = c("panel_id" = "gene_panel")) %>%
    select(-"entrez_ids_in_panel")

  user_data_genes <- gsub(".fus|.Del|.Amp|.cna", "", colnames(gene_binary))

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
                                   annotate_specific_panel,
                                   gene_binary = gene_binary)

  return(annotated_data)
}


#' Utility function  to insert NA's According to Panel
#'
#' @param gene_binary a processed binary matrix
#' @param panel_id name of gene panel
#' @param samples_in_panel samples to be annotated for each panel
#' @param na_genes genes to make NA
#' @param ... other args passed
#'
#' @return an annotated data frame
#' @export
#'
#'
annotate_specific_panel <- function(gene_binary,
                           panel_id,
                           samples_in_panel,
                           na_genes, ...) {

  mut_sub <- gene_binary[samples_in_panel, ]
  mut_sub[,stats::na.omit(match(na_genes, colnames(mut_sub)))] <- NA

  return(mut_sub)

}





