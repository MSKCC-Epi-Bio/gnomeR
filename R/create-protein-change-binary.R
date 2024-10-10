# Protein Change Binary Matrix Function --------------------------

#' Create a binary matrix of protein changes in samples
#'
#' `create_protein_change_binary()` enables creation of a binary matrix of
#' individual protein changes from a mutation file with a predefined list of
#' samples (rows are samples and columns are protein changes).
#'
#' @param mutation A data frame of mutation information in MAF format. Required columns
#' are sample IDs, genes, and protein change names.
#' @param samples A character vector specifying the samples to be included in the output matrix.
#' Default is NULL in which case the function will return all the samples found in the mutation file.
#' If you specify a vector of samples that contain samples not in the mutation data frame, the samples
#' missing from the data file will be removed from the result.
#' @param gene A character vector specifying the target gene(s). Default is NULL and
#' protein changes found in any genes in the mutation file will be returned.
#' @param protein.change A character vector specifying the protein change(s).
#' Default is NULL and function will return all protein changes in the genes.
#'
#' @return A binary protein change data frame (rows are samples and columns are
#' protein changes). The columns names are in the form of `<hugoSymbol>_<proteinChange>`.
#'
#' Values in the data frame are binary indicators `1` (present) and `0` (absent).
#'
#' @export
#' @examples
#' library(gnomeR)
#' all <- create_protein_change_binary(mutation = gnomeR::mutations)
#'
#' gene <- c('EGFR','TP53')
#'
#' double.mut <- create_gene_binary(
#'  mutation = gnomeR::mutations, gene = gene
#' )
#'
#' @import dplyr
create_protein_change_binary <- function(mutation,
                                         samples=NULL,
                                         gene=NULL,
                                         protein.change=NULL) {

  # Check arguments -----------------------------

  ## mutation -------
  ## Check that mutation input is not empty
  if (is.null(mutation)) {
    cli::cli_abort("You must provide the following argument: {.code mutation}.")
  }
  ## Check that mutation input is data.frame
  is_df <- 'data.frame' %in% class(mutation)
  if (!is_df) {
    cli::cli_abort("{.code {not_df}} must be a data.frame")
  }

  ##  samples ------
  if (!(is.null(samples) | is.character(samples))) {
    cli::cli_abort("{.code samples} must be a character vector or `NULL`")
  }

  ##  gene ------
  if (!(is.null(gene) | is.character(gene))) {
    cli::cli_abort("{.code gene} must be a character vector or `NULL`")
  }

  ##  protein changes ------
  if (!(is.null(protein.change) | is.character(protein.change))) {
    cli::cli_abort("{.code samples} must be a character vector or `NULL`")
  }


  # Clean up data input --------------------------

  ## Mutation dataframe -------------------
  ## standardize columns names
  mutation <- switch(!is.null(mutation),
                     .clean_and_check_cols(
                       df_to_check = mutation,
                       required_cols = c("sample_id", "hugo_symbol")))

  ## Protein change symbol ----------------
  ## remove the "p." symbol (and others) at the beginning if exist
  protein.change <- protein.change %||% mutation$hgv_sp_short %>%
    as.character() %>%
    unique()

  protein.change <- as.vector(sapply(protein.change,
                                     function(x) stringr::str_replace(x,
                                                                      "\\bp.|\\bc.|\\bg.|\\bm.|\\br.",
                                                                      ""),
                                     simplify = T))

  ## Final Gene List ----------------------
  gene <- gene %||% mutation$hugo_symbol %>%
    as.character() %>%
    unique()


  ## Final Sample List --------------------

  # get whole list of samples in the mutation data frame
  samples_in_data <- mutation$sample_id %>% as.character() %>% unique()
  samples = unique(samples)

  # if provided samples not found in the mutation dataset
  if (!is.null(samples) & all(!(samples %in% samples_in_data))) {

    cli::cli_abort("None of your selected {.code samples} have alterations in your data.")

  } else if( length(setdiff(samples, samples_in_data))!=0 ) {
    null_samples = setdiff(samples, samples_in_data)
    cli::cli_alert_warning("Samples {null_samples} not found in the {.field mutation} data frame. We will remove them from the resulting data frame.")

    samples_final = intersect(samples, samples_in_data)

  } else if (is.null(samples)) {

    # if samples not passed we will infer it from data frame
    cli::cli_alert_warning("{.code samples} argument is {.code NULL}. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in the {.field mutation} data frame.")

    # If user doesn't pass a vector, use samples in files as final sample list
    samples_final <- samples %||%
      samples_in_data %>%
      unique()

  } else {
    samples_final = samples
  }


  # Create mutation binary matrix -------------------------

  protein.change.binary <- mutation %>%
    ## filter for the target gene and protein change
    filter(hugo_symbol  %in% gene &
           hgv_sp_short %in% protein.change) %>%
    mutate(point_mutation = paste(hugo_symbol, hgv_sp_short, sep="_"))  %>%

    select(sample_id,
           hugo_symbol,
           point_mutation)  %>%

    ## select unique combination of patient and point mutation (optional)
    group_by(sample_id, hugo_symbol, point_mutation)  %>%
    filter(row_number()==1)  %>%
    ungroup() %>%
    mutate(fl = 1) %>%

    ## convert into wide format with binary indicators for mutations
    tidyr::pivot_wider(id_cols = "sample_id", names_from = "point_mutation", values_from  = "fl",
                       values_fill = 0,
                       # names_glue = rlang::eval_tidy(names_glue)
    )  %>%
    right_join(
      data.frame(sample_id = samples_final),
      by='sample_id'
    ) %>%
    replace(is.na(.), 0)

  return(protein.change.binary)
}
