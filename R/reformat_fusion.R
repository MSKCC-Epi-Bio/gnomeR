#' Enables users to reformat fusions files so that each fusion is listed as one row with two hugo-symbol
#' sites instead of two rows, one for each site. This is the required format for the `create_gene_binary` function.
#'
#' @param fusions a data frame of fusion/structural variants that occur in a cohort. There should be a `sample_id`, `hugo_symbol`,
#' and `fusion` column at minimum. Intragenic/intergenic fusions will have one row. Any two gene fusions will
#' have two rows. See `gnomeR::sv_long` for an example.
#'
#' @return a data frame with `sample_id`, `site1hugo_symbol`, and `site2hugo_symbol` and `fusion` columns. This should match the format
#' of the `gnomeR::sv` dataset.
#' @export
#' @examples
#'
#' sv_long1 <- gnomeR::sv_long %>%
#'   rename_columns() %>%
#'   reformat_fusion()
#'
#' head(sv_long1)
#'
#' @import dplyr
#' @import stringr
#' @import tidyr



# Function ----------------------------------------------------------------

reformat_fusion <- function(fusions) {

  # Checks --------------------------------------------------------------------
  if (!is.data.frame(fusions)) {
    cli::cli_abort("{.code fusion} must be a data.frame")
  }

  fusions <- gnomeR::rename_columns(fusions)
  .check_required_cols(fusions, c("fusion", "sample_id", "hugo_symbol"))

  # Clean Strings -------------------------------------------------------------
  fusions_sep <- fusions %>%
    mutate(
      # remove leading space in fusion var
      fusion = stringr::str_trim(.data$fusion),
      # ensure there are no spaces between hyphen and next letter (problem with - Archer)
      fusion = stringr::str_replace_all(.data$fusion, "- ", "-"),
      # remove "Archer"
      fusion = stringr::str_remove_all(.data$fusion, "-Archer"))


  # Replace Known Hugo Exceptions -------------------------------------------

  # check for cells with more than one hyphen
  any_double_dash <- (stringr::str_count(fusions_sep$hugo_symbol, "-") == 1 |
                        stringr::str_count(fusions_sep$fusion, "-") > 1)

  # list of known genes with a hyphen to replace with underscore in data
  exceptions <- c(
    # generalized edits
    "-AS" = "_AS",
    "HLA-" = "HLA_",
    "NKX2-" = "NKX2_",
    "NKX3-" = "NKX3_",

    # specific edits
    "SOX2-OT"= "SOXOT",
    "MIR365-2"= "MIR365_2",
    "LINC-PINT"= "LINC_PINT",
    "NKX2-8"= "NKX2_8",
    "PMF1-BGLAP"= "PMF1_BGLAP",
    "BIVM-ERCC5"= "BIVM_ERCC5",
    "PRH1-PRR4"= "PRH1_PRR4",
    "RNU6-19P"= "RNU6_19P",
    "BRWD1-IT2"= "BRWD1_IT2",
    "CTD-2151A2.1"= "CTD_2151A2.1")

  fusions_sep$fusion[any_double_dash] <- stringr::str_replace_all(
    fusions_sep$fusion[any_double_dash], exceptions)

  # find any remaining hugo_symbols that have 2 hyphens
  any_gr3 <- fusions_sep$fusion[stringr::str_count(fusions_sep$fusion, "-") > 1]

  if(length(any_gr3 > 0)) {

    cli::cli_abort(c("Unable to process the following fusions in the 'fusion'",
                       "col: {any_gr3}. Please reformat these as `HUGO1-HUGO2`. Hugo names cannot have '-' only '_'"
      ))
    }



  # Split Into Hugo Symbol 1 & 2  ---------------------------------------------

  # only retain string before and after first dash, then remove intergenic/intragenic
  fusions_sep <- fusions_sep %>%
    mutate(fusion_orig = .data$fusion) %>%
    mutate(
      fusion = stringr::str_extract(.data$fusion, "[[:alnum:]_]+-[[:alnum:]_]+")) %>%
    mutate(
      fusion = stringr::str_remove_all(.data$fusion,
                              "-INTERGENIC|-intergenic|-INTRAGENIC|-intragenic"))

  # split into hugo symbol columns
  fusions_sep <- fusions_sep %>%
    separate_wider_delim("fusion",
                         names = paste0("site", 1:4, "hugo_symbol"),
             delim = "-",
             cols_remove = FALSE,
             too_many = "drop",
             too_few = "align_start")



  # filters out any mistakes because this should never be in the first site position
  fusions_sep <- fusions_sep %>%
    select(
      "sample_id", "hugo_symbol", "site1hugo_symbol", "site2hugo_symbol", "fusion"
    ) %>%
    filter(!(.data$site1hugo_symbol %in% c("repeat", "insufficient")))


  # There are cases where site 1 and 2 were flipped for a sample_id and listed x2
  # ex: TP53-APC vs APC-TP53
  # Here we sort each hugo alphabetically, then paste
  fusions_sep <- fusions_sep %>%
    mutate(fusions_ordered = purrr::pmap_chr(list(.data$site1hugo_symbol, .data$site2hugo_symbol),
                                             ~paste0(sort(c(.x, .y)), collapse = "-")))

  fusions_unq <- fusions_sep %>%
    select("sample_id", "fusions_ordered") %>%
    distinct()

  fusions_unq <- fusions_unq %>%
    separate_wider_delim("fusions_ordered",
                         names = paste0("site_", 1:4, "_hugo_symbol"),
                         delim = "-",
                         cols_remove = FALSE,
                         too_many = "drop",
                         too_few = "align_start") %>%
    rename("fusion" = "fusions_ordered")

  fusions_unq
}


