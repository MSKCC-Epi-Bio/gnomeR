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

  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
  if(!is.character(maf$Hugo_Symbol)) maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
  if(!is.character(maf$Tumor_Sample_Barcode)) maf$Tumor_Sample_Barcode <-
      as.character(maf$Tumor_Sample_Barcode)

  # get table of gene aliases
  alias_table <- tidyr::unnest(impact_genes, cols = alias) %>%
    select(hugo_symbol, alias)

  # recode aliases
  maf$Hugo_Symbol_Old <- maf$Hugo_Symbol
  maf$Hugo_Symbol <- purrr::map_chr(maf$Hugo_Symbol, ~resolve_alias(.x,
                                      alias_table = alias_table))

  message <- maf %>%
    filter(Hugo_Symbol_Old != Hugo_Symbol) %>%
    select(Hugo_Symbol_Old, Hugo_Symbol) %>%
    distinct()

  if(nrow(message) > 0) {
    warning(paste0("To ensure gene with multiple names/aliases are correctly grouped together, the
    following genes in your maf dataframe have been recoded: \n",
                   purrr::map2(message$Hugo_Symbol_Old,
                        message$Hugo_Symbol,
                        ~paste0(.x, " recoded to ", .y, " \n"))))
  }

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
#' @param gene_to_check
#' @param alias_table
#'
#' @return if the accepted hugo symbol is input, it is returned back.
#' If an alias name is provided, the more common name/more up to date name is returned
#' @export
#'
#' @examples
#' resolve_alias("KMT2D", alias_table = tidyr::unnest(impact_genes, cols = alias))
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
