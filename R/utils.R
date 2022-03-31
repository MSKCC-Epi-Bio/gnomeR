#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param ... Further arguments parsed through binmat() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#' check_maf_input(mut,recode.aliases = TRUE)
#'
check_maf_input <- function(maf, ...)  {

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



  # * Recode Gene Aliases---------------------

  if(arguments$recode.aliases == TRUE) {

  # get table of gene aliases (internal data)
    alias_table <- tidyr::unnest(impact_gene_info, cols = .data$alias) %>%
      dplyr::select(.data$hugo_symbol, .data$alias)

    # recode aliases
    maf$Hugo_Symbol_Old <- maf$Hugo_Symbol
    maf$Hugo_Symbol <- purrr::map_chr(maf$Hugo_Symbol, ~resolve_alias(gene_to_check = .x,
                                                                      alias_table = alias_table))

    message <- maf %>%
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
      following genes in your maf dataframe have been recoded (you can supress this with {.code recode.aliases = FALSE}):",
      vec_recode))

    }
  }

  return(maf)
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
