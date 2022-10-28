#' tbl_genomic
#'
#' This function will select genes based on user inputs or on frequency counts and then
#' will pass the data.frame to `gtsummary::tbl_summary()`. You can specify a `by` variable and other
#' parameters that are accepted by `gtsummary::tbl_summary()`. Note the `by` variable must be merged on to
#' onto the data before using the `by` parameter in the function.
#'
#' @param gene_binary data.frame of genetic samples
#' @param freq_cutoff A number 0 to 1 representing the minimum percent frequency you are using to select gene's for analysis.
#' Frequencies can be calculated at gene level, or alteration level (see `freq_cutoff_by_gene`).
#' @param freq_cutoff_by_gene Logical indicating whether gene selection based on frequency % should
#' be calculated at the gene level, or the alteration level. Default is TRUE, indicating all
#' alterations (e.g. TP53/TP53.Del/TP53.Del/TP53.fus) will be aggregated at gene level to
#' determine cutoff frequencies and all alteration types for top genes will be returned.
#' @param gene_subset Specific genes you want to summarize. This takes precedent over the `freq_cutoff` parameter and will
#' return all alterations associated with that gene (including mutations, CNA, fusions).
#' @param by A variable to be passed to `gtsummary::tbl_summary()`'s by parameter
#' @param ... Additional parameters that can be passed to `gtsummary::tbl_summary()`. To access the additional parameters you need to load `gtsummary`.
#'
#' @return A `tbl_summary()` object
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' samples <- unique(mutations$sampleId)[1:10]
#'
#' gene_binary <- create_gene_binary(
#'   samples = samples,
#'   mutation = gnomeR::mutations,
#'   cna = gnomeR::cna,
#'   mut_type = "somatic_only", snp_only = FALSE,
#'   specify_panel = "no"
#' )
#' tb1 <- tbl_genomic(gene_binary = gene_binary, freq_cutoff = .05)
#' tb2 <- tbl_genomic(gene_binary = gene_binary, gene_subset = c("BRCA1", "FGFR2"))
#'
#' gene_binary <- gene_binary %>%
#'   mutate(sex = sample(
#'     x = c("M", "F"),
#'     size = nrow(gene_binary), replace = TRUE
#'   ))
#' library(gtsummary) # Need to load gtsummary to access additional arguments in `...`
#'
#' t1 <- tbl_genomic(
#'   gene_binary = gene_binary,
#'   by = sex,
#'   freq_cutoff = .2,
#'   freq_cutoff_by_gene = FALSE,
#'   statistic = list(all_categorical() ~"{n}")
#' )
#'
tbl_genomic <- function(gene_binary,
                        freq_cutoff = NULL,
                        freq_cutoff_by_gene = TRUE,
                        gene_subset = NULL,
                        by = NULL, ...) {

  # check arguments & prep data ------------------------------------------------

  if (!inherits(gene_binary, "data.frame")) {
    stop("`gene_binary=` argument must be a tibble or data frame.", call. = FALSE)
  }

  if (!("sample_id" %in% names(gene_binary))) {
    gene_binary <- tibble::rownames_to_column(gene_binary, var = "sample_id")
  }

  if (!is.null(freq_cutoff) && (freq_cutoff < 0 || freq_cutoff > 1)) {
    cli::cli_abort("Please select a {.code freq_cutoff} value between {.code 0} and {.code 1}")
  }

  if (is.null(freq_cutoff) & is.null(gene_subset)) {
    cli::cli_alert("Please note neither {.code freq_cutoff} or {.code gene_subset} had inputs. By default {.code freq_cutoff} will be set to 0.1")
    freq_cutoff <- .1
  }

  by <-
    .select_to_varnames({{ by }},
      data = gene_binary,
      arg_name = "by", select_single = TRUE
    )


  # check & assign gene subset -------------------------------------------------

  # if user passes gene_subset, we will add sufix
  gene_subset <- gene_subset %>%
    purrr::when(
      is.null(.) ~ .,
      !is.character(gene_subset) ~
        cli::cli_abort("Please supply a character vector for {.code gene_subset}"),
      length(.[(. %in% colnames(gene_binary))]) == 0 ~
        cli::cli_abort("No genes specified in {.code gene_subset} are in your gene_binary"),
      str_detect(., ".Amp|.Del|.fus|.cna") ~
        cli::cli_abort(
          "Detected one of the following in {.code gene_subset}: {.code '.Amp|.Del|.fus|.cna'} You may
          only pass gene names (eg. 'TP53'). To only include specific alterations, consider {.code dplyr::select(df, <alterations>)}
          before passing to {.code tbl_genomic()}"
        ),
      freq_cutoff > 0 ~ {
        cli::cli_inform("You've supplied both {.code gene_subset} and {.code freq_cutoff}.
                       {.code freq_cutoff} parameter will be ignored")
        return(.)
      },

      # return only genes found in your data
      length(.[!(. %in% colnames(gene_binary))]) > 0 ~ {
        cli::cli_warn("The following of {.code gene_subset} are not in your data: {.code { .[!(. %in% colnames(gene_binary))]}}")
        return(.[(. %in% colnames(gene_binary))])
      },
      TRUE ~ .
    )

  gene_subset <- switch(!is.null(gene_subset),
    c(
      gene_subset,
      paste0(gene_subset, ".Amp"),
      paste0(gene_subset, ".Del"),
      paste0(gene_subset, ".fus"),
      paste0(gene_subset, ".cna")
    )
  )

  # Calc Gene Frequencies (if gene_subset is NULL) --------------------------
  gene_subset <- gene_subset %||% {
    gene_binary %>%
      select(-all_of(by)) %>%
      # if freq should be calc at gene level- simplify matrix first
      # todo- if simplify matrix already called- avoid this!
      purrr::when(
        freq_cutoff_by_gene ~ summarize_by_gene(.),
        TRUE ~ .
      ) %>%
      ungroup() %>%
      tidyr::pivot_longer(-.data$sample_id) %>%
      distinct() %>%
      group_by(.data$name) %>%
      summarise(
        sum = sum(.data$value, na.rm = TRUE),
        count = nrow(gene_binary) - sum(is.na(.data$value)),
        num_na = sum(is.na(.data$value))
      ) %>%
      mutate(perc = .data$sum / .data$count) %>%
      filter(.data$perc >= freq_cutoff) %>%
      arrange(desc(.data$perc)) %>%
      pull("name")
  }

  if (length(gene_subset) < 1) {
    cli::cli_abort("No genes in data set match your filter criteria (see {.code freq_cutoff})")
  }

  # if freq cutoff by gene
  gene_subset <- gene_subset %>%
    purrr::when(
      freq_cutoff_by_gene ~ c(
        .,
        paste0(., ".Amp"),
        paste0(., ".Del"),
        paste0(., ".fus"),
        paste0(., ".cna")
      ),
      TRUE ~ .
    )


  # Select Genes and Make Table-----------------------------------------------

  table_data <- gene_binary %>%
    select(all_of(by), any_of(c(gene_subset)))



   table_data %>%
    gtsummary::tbl_summary(by = by,...)

  # should we split results by MUT/CNA/Fusion? if not at least have to fix arranging of  these
  # also its confusing when freq_by_gene is TRUE and you see low freq alts in table

}
