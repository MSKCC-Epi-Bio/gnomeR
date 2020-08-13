#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param maf
#'
#' @return a corrected maf file or an error if problems with maf
#' @export
#'
#' @examples
#'
#'check_maf_input(mut)
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
  if(!is.character(maf$Tumor_Sample_Barcode)) maf$Tumor_Sample_Barcode <- as.character(maf$Tumor_Sample_Barcode)
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

  return(maf)
}


calc_variant_prop <- function(maf) {
    maf %>%
    group_by(.data$Tumor_Sample_Barcode) %>%
    mutate(totalMut = n()) %>%
    ungroup() %>%
    group_by(.data$Tumor_Sample_Barcode,.data$Variant_Classification) %>%
    summarise(N=n(),
              varProp = .data$N/unique(.data$totalMut))

}

  # function to extract SNV class from character string cols
  substrRight <- function(x, n) {
    x <- as.character(x)
    substr(x, nchar(x) - n + 1, nchar(x))
  }
