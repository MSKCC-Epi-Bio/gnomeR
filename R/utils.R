#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param maf Raw maf dataframe containing alteration data
#'
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#'
#' check_maf_input(mut)
#'
check_maf_input <- function(maf)  {

  # data checks for maf files
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

    if(is.na(match("Reference_Allele",colnames(maf))) || is.na(match("Tumor_Seq_Allele2",colnames(maf)))){
      warning("The MAF file inputted is missing a mutation status column (Variant_Type).
            The variant type couldn't be inferred from the data (columns Reference_Allele and Tumor_Seq_Allele2 are required)
              and it will be assumed that all variants are SNPs.")
      maf$Variant_Type <- rep("SNPs",nrow(maf))
    }

    else{
      warning("The MAF file inputted is missing a mutation status column (Variant_Type).
            The variant type will be inferred (we recommend that the user fixes this).")
      maf <- maf %>%
        mutate(
          Reference_Allele = as.character(Reference_Allele),
          Tumor_Seq_Allele2 = as.character(Tumor_Seq_Allele2),
          Variant_Type = case_when(
            Reference_Allele %in% c("A","T","C","G") &
              Tumor_Seq_Allele2 %in% c("A","T","C","G") ~ "SNP",
            nchar(Tumor_Seq_Allele2) < nchar(Reference_Allele) |
              Tumor_Seq_Allele2 == "-" ~ "DEL",
            Reference_Allele == "-" |
              nchar(Tumor_Seq_Allele2) > nchar(Reference_Allele) ~ "INS",
            nchar(Reference_Allele) == 2 & nchar(Tumor_Seq_Allele2) == 2 ~ "DNP",
            nchar(Reference_Allele) == 3 & nchar(Tumor_Seq_Allele2) == 3 ~ "TNP",
            nchar(Reference_Allele) > 3 & nchar(Tumor_Seq_Allele2) == nchar(Reference_Allele) ~ "ONP",
            TRUE ~ "Undefined"
          ))
    }
  }

  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
  if(!is.character(maf$Hugo_Symbol)) maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
  if(!is.character(maf$Tumor_Sample_Barcode)) maf$Tumor_Sample_Barcode <-
      as.character(maf$Tumor_Sample_Barcode)

  # get table of gene aliases
  alias_table <- tidyr::unnest(impact_gene_info, cols = alias) %>%
    dplyr::select(hugo_symbol, alias)

  # recode aliases
  maf$Hugo_Symbol_Old <- maf$Hugo_Symbol
  maf$Hugo_Symbol <- purrr::map_chr(maf$Hugo_Symbol, ~resolve_alias(.x,
                                                                    alias_table = alias_table))

  message <- maf %>%
    dplyr::filter(Hugo_Symbol_Old != Hugo_Symbol) %>%
    dplyr::select(Hugo_Symbol_Old, Hugo_Symbol) %>%
    dplyr::distinct()

  if(nrow(message) > 0) {
    warning(paste0("To ensure gene with multiple names/aliases are correctly grouped together, the
    following genes in your maf dataframe have been recoded: \n",
                   purrr::map2(message$Hugo_Symbol_Old,
                               message$Hugo_Symbol,
                               ~paste0(.x, " recoded to ", .y, " \n"))))
  }

  # maf <- maf %>%
  #   mutate(Hugo_Symbol = gsub("-",".",as.character(Hugo_Symbol)))
  return(maf)
}


#' Check Columns of MAF
#'
#' @param maf Raw maf dataframe containing alteration data
#' @param col_to_check Which column to check
#'
#' @return
#' @export
#'
#' @examples
#' check_maf_column(mut, "Variant_Classification")
#'
check_maf_column <- function(maf, col_to_check) {
  if(!(col_to_check %in% colnames(maf))) {
    stop(paste("The MAF file inputted is missing the following column:", col_to_check))
  }
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
      filter(alias == gene_to_check) %>%
      pull(hugo_symbol) %>%
      first()

  } else {
    gene_to_check
  }
}
