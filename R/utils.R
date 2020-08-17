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
