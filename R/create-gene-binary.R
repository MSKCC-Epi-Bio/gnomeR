# Main Binary Matrix Function ------------------

#' Enables creation of a binary matrix from a mutation, fusion or CNA file with
#' a predefined list of samples (rows are samples and columns are genes)
#' @param samples a character vector specifying which samples should be included in the resulting data frame.
#' Default is NULL is which case all samples with an alteration in the mutation, cna or fusions file will be used. If you specify
#' a vector of samples that contain samples not in any of the passed genomic data frames, 0's (or NAs when appropriate if specifying a panel) will be
#' returned for every column of that patient row.
#' @param mutation A data frame of mutations in the format of a maf file.
#' @param mut_type The mutation type to be used. Options are "omit_germline", "somatic_only", "germline_only" or "all". Note "all" will
#' keep all mutations regardless of status (not recommended). Default is omit_germline which includes all somatic mutations, as well as any unknown mutation types (most of which are usually somatic)
#' @param snp_only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include_silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @param fusion A data frame of fusions. If not NULL the outcome will be added to the matrix with columns ending in ".fus".
#' Default is NULL.
#' @param cna A data frame of copy number alterations. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
#' Default is NULL.
#' @param high_level_cna_only If TRUE, only deep deletions (-2, -1.5) or high level amplifications (2) will be counted as events
#' in the binary matrix. Gains (1) and losses (1) will be ignored. Default is `FALSE` where all CNA events are counted.
#' @param specify_panel Default is `"no"` where no panel annotation is done. Otherwise pass a character vector of length 1 with a
#' panel id (see `gnomeR::gene_panels` for available panels), or `"impact"` for automated IMPACT annotation.
#' Alternatively, you may pass a data frame of `sample_id`-`panel_id` pairs specifying panels for each sample for
#' which to insert NAs indicating genes not tested. See below for details.
#' @param recode_aliases Default is `"impact"` where function will check for IMPACT genes that may go by more than 1 name in your data and replace the alias name with the standardized gene name (see `gnomeR::impact_alias_table` for reference list).
#' If `"no"`, no alias annotation will be performed. If `"genie"`, an alias table with GENIE BPC genes will be used to check (see `gnomeR::genie_alias_table` for reference list).
#' Alternatively, you may pass a custom alias list as a data frame with columns `hugo_symbol` and `alias` specifying a custom alias table to use for checks. See below for details.
#'
#'
#' @section `specify_panel` argument:
#'    - If `specify_panel = "no"` is passed (default) data will be returned as is without any additional NA annotations.
#'    - If a single panel id is passed (e.g. `specify_panel = "IMPACT468"`), all genes in your data that are not tested on that panel will be set to
#' `NA` in results for all samples (see gnomeR::gene_panels to see which genes are on each supported panels).
#'    - If `specify_panel = "impact"` is passed, impact panel version will be inferred based on each sample_id (based on `IMX` nomenclature) and NA's will be
#' annotated accordingly for each sample/panel pair.
#'    - If you wish to specify different panels for each sample, pass a data frame (with all samples included) with columns: `sample_id`, and `panel_id`. Each sample will be
#' annotated with NAs according to that specific panel. If a sample in your data is missing from the `sample_id` column in the
#' `specify_panel` dataframe, it will be returned with no annotation (equivalent of setting it to "no").
#'
#' @section `recode_aliases` argument:
#'    - If `recode_aliases = "impact"` is passed (default), function will use `gnomeR::impact_alias_table` to find and replace any non-standard hugo symbol names with their
#'    more common (or more recent) accepted gene name.
#'    - If `recode_aliases = "genie"` is passed, function will use `gnomeR::genie_alias_table` to find and replace any non-standard hugo symbol names with their
#'    more common (or more recent) accepted gene name.
#'    - If `recode_aliases = "no"` is passed, data will be returned as is without any alias replacements.
#'    - If you have a custom table of vetted aliases you wish to use, you can pass a data frame with columns: `hugo_symbol`, and `alias`.
#'      Each row should have one gene in the `hugo_symbol` column indicating the accepted gene name, and one gene in the `alias` column indicating an alias
#'      you want to check for and replace. If a gene has multiple aliases to check for, each should be represented in its own separate row.
#'      See `gnomeR::impact_alias_table` for an example of accepted data formatting.
#'
#'
#' @return a data frame with sample_id and alteration binary columns with values of 0/1
#' @export
#' @examples
#' mut.only <- create_gene_binary(mutation = gnomeR::mutations)
#'
#' samples <- gnomeR::mutations$sampleId
#'
#' bin.mut <- create_gene_binary(
#'   samples = samples, mutation = gnomeR::mutations,
#'   mut_type = "omit_germline", snp_only = FALSE,
#'   include_silent = FALSE
#' )
#'
#' @import dplyr
#' @import stringr

create_gene_binary <- function(samples = NULL,
                               mutation = NULL,
                               mut_type = c("omit_germline", "somatic_only", "germline_only", "all"),
                               snp_only = FALSE,
                               include_silent = FALSE,
                               fusion = NULL,
                               cna = NULL,
                               high_level_cna_only = FALSE,
                               specify_panel = "no",
                               recode_aliases = "impact") {
  pathways <- gnomeR::pathways
  gene_panels <- gnomeR::gene_panels

  # Check Arguments ------------------------------------------------------------

  if (is.null(mutation) && is.null(fusion) && is.null(cna)) {
    cli::cli_abort("You must provide at least one of the three following arguments: {.code mutation}, {.code fusion} or {.code cna}.")
  }

  # Check that mutation, fusion, cna is data.frame
  is_df <- purrr::map(
    list(mutation = mutation, fusion = fusion, cna = cna),
    ~ dplyr::case_when(
      !is.null(.x) ~ "data.frame" %in% class(.x)
    )
  ) %>%
    purrr::compact()

  not_df <- names(is_df[which(is_df == FALSE)])

  if (length(not_df) > 0) {
    cli::cli_abort("{.code {not_df}} must be a data.frame")
  }

  # * mut_type-----
  mut_type <- match.arg(mut_type)

  # * Specify Panel --------
  # must be a known character or data frame with specified column

  # make tibbles into data.frames - idk if this is needed, could change switch to ifelse I think a alternative
  if ("tbl" %in% class(specify_panel)) {
    specify_panel <- as.data.frame(specify_panel)
  }

  specify_panel <-
    switch(class(specify_panel),
      "character" = {
        choices_arg <- c("no", "impact", "IMPACT", gene_panels$gene_panel)
        match.arg(specify_panel, choices = choices_arg)
      },
      "data.frame" = {
        # check for correct column names
        if (!("sample_id" %in% names(specify_panel)) | !("panel_id" %in% names(specify_panel))) {
          cli::cli_abort(c(
            "Dataframe passed to {.var specify_panel} must have columns for ",
            "{.code sample_id} and {.code panel_id}."
          ))
        }

        if (any(is.na(specify_panel$panel_id))) {
          cli::cli_abort("Some {.field panel_id} values in {.code sample_panel_pair} df are {.code NA}. Please explicitely indicate {.code no} for those samples instead if you wish to skip annotating these.")
        }

        if (length(setdiff(c(specify_panel$panel_id), c(gene_panels$gene_panel, "no"))) > 0) {
          cli::cli_abort("Panels not known: {.val {setdiff(c(specify_panel$panel_id), c(gene_panels$gene_panel, 'no'))}}. See {.code  gnomeR::gene_panels} for known panels, or skip annotation with {.code specify_panel = 'no'} or indicating {.code 'no'} for those samples in {.field panel_id} column of sample_id-panel_id pair data frame")
        }
        specify_panel
      },
      cli::cli_abort("{.code specify_panel} must be a character vector of length 1 or a data frame.")
    )


  # * Mutation  checks  --------

  # standardize columns names
  mutation <- switch(!is.null(mutation),
    sanitize_mutation_input(
      mutation = mutation,
      include_silent = include_silent
    )
  )
  names_mut_dict <- attr(mutation, "names_mut_dict")

  # * Fusion checks  ----------
  fusion <- switch(!is.null(fusion),
    sanitize_fusion_input(fusion)
  )

  # * CNA checks  ------------
  cna <- switch(!is.null(cna),
    {
      sanitize_cna_input(cna)
    }
  )


  #  Make Final Sample List ----------------------------------------------------

  samples_in_data <-
    c(mutation$sample_id, fusion$sample_id, cna$sample_id) %>%
    as.character() %>%
    unique()


  if(!is.null(samples) & all(!(samples %in% samples_in_data))) {
    cli::cli_abort("None of your selected {.code samples} have alterations in your data. ")
  }

  # if samples not passed we will infer it from data frames
  samples %||%
    cli::cli_alert_warning("{.code samples} argument is {.code NULL}. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in {.field mutation}, {.field fusion} or {.field cna} data frames")

  # If user doesn't pass a vector, use samples in files as final sample list
  samples_final <- samples %||%
    samples_in_data

  # Recode Aliases -----------------------------------------------------------

  # Fusions - create long version with event split by two involved genes
  if(!is.null(fusion)) {
    fusion <- fusion %>%
      select(
        "sample_id",
        "site_1_hugo_symbol",
        "site_2_hugo_symbol"
      ) %>%
      tidyr::pivot_longer(-"sample_id", values_to = "hugo_symbol") %>%
      select("sample_id", "hugo_symbol")
  }

  if (recode_aliases != "no") {

    all_alias_warnings <- c()

    if(!is.null(mutation)) {
      q_mut <- recode_alias(mutation,
                            alias_table = recode_aliases, supress_warnings = TRUE)
      mutation <- q_mut$genomic_df
      q_mut_warn <- q_mut$aliases_in_data
      all_alias_warnings <- c(all_alias_warnings, q_mut_warn)
    }

    if(!is.null(cna)) {
      q_cna <- recode_alias(cna, alias_table = recode_aliases, supress_warnings = TRUE)
      cna <- q_cna$genomic_df
      q_cna_warn <- q_cna$aliases_in_data
      all_alias_warnings <- c(all_alias_warnings, q_cna_warn)
    }

    if(!is.null(fusion)) {
      q_fus <- recode_alias(fusion, alias_table = recode_aliases, supress_warnings = TRUE)
      fusion <- q_fus$genomic_df
      q_fus_warn <- q_fus$aliases_in_data
      all_alias_warnings <- c(all_alias_warnings, q_fus_warn)
    }

    all_alias_warnings <- unique(all_alias_warnings)

    if (length(all_alias_warnings) > 0) {
      cli::cli_warn(c(
        "To ensure gene with multiple names/aliases are correctly grouped together, the
        following genes in your dataframe have been recoded (if you are running {.code create_gene_binary()}
        you can prevent this with {.code alias_table = FALSE}):",
        all_alias_warnings
      ))
    }
  }


  # Binary matrix for each data type ------------------------------------------

  # create quiet versions to catch and combine messages
  mutation_binary_df <- switch(!is.null(mutation),
    .mutations_gene_binary(
      mutation = mutation,
      samples = samples_final,
      mut_type = mut_type,
      snp_only = snp_only,
      include_silent = include_silent,
      specify_panel = specify_panel
    )
  )

  # fusions
  fusion_binary_df <- switch(!is.null(fusion),
    .fusions_gene_binary(
      fusion = fusion,
      samples = samples_final,
      specify_panel = specify_panel
    )
  )

  # cna
  cna_binary_df <- switch(!is.null(cna),
    .cna_gene_binary(
      cna = cna,
      samples = samples_final,
      specify_panel = specify_panel,
      high_level_cna_only = high_level_cna_only
    )
  )

  # put them all together
  df_list <- list(mutation_binary_df, fusion_binary_df, cna_binary_df)

  all_binary <- purrr::reduce(df_list[!sapply(df_list, is.null)], # remove null if present
    full_join,
    by = "sample_id"
  ) %>%
    mutate(across(setdiff(everything(), "sample_id"), .fns = function(x) {
      ifelse(is.na(x), 0, x)
    }))

  # add in any samples with no mutations
  if (!is.null(samples)) {
    no_alt_samples <- setdiff(samples_final, all_binary$sample_id)

    if (length(no_alt_samples) > 0) {
      add_no_alt_samples <-
        data.frame(matrix(0, ncol = ncol(all_binary), nrow = length(no_alt_samples)))

      names(add_no_alt_samples) <- names(all_binary)
      add_no_alt_samples$sample_id <- no_alt_samples

      all_binary <- bind_rows(all_binary, add_no_alt_samples)
      all_binary <- all_binary[match(samples_final, all_binary$sample_id), ]

    }
  }

  # Platform-specific NA Annotation ------

  # we've already checked the arg is valid
  # If character, make into data frame sample-panel pair to input in function
  if (is.character(specify_panel)) {
    sample_panel_pair <- switch(specify_panel,
      "impact" = specify_impact_panels(all_binary),
      "no" = {
        all_binary["sample_id"] %>%
          mutate(panel_id = "no")
      },
      all_binary["sample_id"] %>%
        mutate(panel_id = specify_panel)
    )
    # create data frame of sample IDs
  } else {
    specify_panel <- specify_panel %>%
      select("sample_id", "panel_id")

    diff_samp <- setdiff(samples_final, specify_panel$sample_id)

    if (length(diff_samp) > 0) {
      # If some samples are not in the specify_panel df, add them as no annotation.
      # TODO Should we add warning?
      add_on <- cbind.data.frame("sample_id" = diff_samp, "panel_id" = rep("no", length(diff)))

      specify_panel <- rbind.data.frame(specify_panel, add_on)
    }

    sample_panel_pair <- specify_panel
  }

  all_binary <- annotate_any_panel(sample_panel_pair, all_binary)


  # Warnings and Attributes --------

  # Throw Message About Empty Columns ------
  all_column_is_na <- names(all_binary)[apply(all_binary, 2, function(x) sum(is.na(x))) == nrow(all_binary)]

  if (length(all_column_is_na) > 0) {
    cli::cli_alert_warning(c(
      "{length(all_column_is_na)} column{?s} {?has/have} all missing values. This may occur when ",
      "there are genes in your data that are not in the specified panels (see `specify_panel` argument)"
    ))
  }

  # return omitted zero  samples as warning/attribute
  samples_no_alts <- setdiff(samples_final, samples_in_data)

  if(length(samples_no_alts) > 0) {
    attr(all_binary, "zero_alteration_samples") <- samples_no_alts

    cli::cli_alert_warning(c("{length(samples_no_alts)} {.code samples} had no alterations ",
    "found in data sets (See {.code attr(<your_df>, 'zero_alteration_samples')} to view). ",
    "These were retained in results as having 0 alterations."))

  }

  return(all_binary)
}


# Mutations Binary Matrix -----------------------------------------------------

#' Make Binary Matrix From Mutation data frame
#'
#' @inheritParams create_gene_binary
#' @keywords internal
#' @return a data frame
#' @export
#'
.mutations_gene_binary <- function(mutation,
                                   samples,
                                   mut_type,
                                   snp_only,
                                   include_silent,
                                   specify_panel) {

  # apply filters --------------

  if (snp_only) {
    mutation <- filter(mutation, .data$variant_type == "SNP")
  }

  if (include_silent == FALSE) {
    mutation <- filter(
      mutation,
      .data$variant_classification != "Silent" |
        is.na(.data$variant_classification)
    )
  }


  switch(mut_type,
    "all" = {
      mutation <- mutation
    },
    "omit_germline" = {
      mutation <- mutation %>%
        filter(.data$mutation_status != "GERMLINE" |
          .data$mutation_status != "germline" | is.na(.data$mutation_status))

      blank_muts <- mutation %>%
        filter(is.na(.data$mutation_status) |
          .data$mutation_status == "" |
          .data$mutation_status == "NA") %>%
        nrow()

      if ((blank_muts > 0)) {
        cli::cli_alert_warning(
          "{(blank_muts)} mutations have {.code NA} or blank in the {.field {names_mut_dict['mutation_status']}} column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.")
      }
    },
    "somatic_only" = {
      mutation <- mutation %>%
        filter(.data$mutation_status == "SOMATIC" |
          .data$mutation_status == "somatic")
    },
    "germline_only" = {
      mutation <- mutation %>% filter(.data$mutation_status == "GERMLINE" |
        .data$mutation_status == "germline")
    }
  )



  mut_bm <- .process_binary(data = mutation, samples = samples, type = "mut")

  return(mut_bm)
}



# Fusions Binary Matrix -----------------------------------------------------

#' Make Binary Matrix From Fusion data frame
#'
#' @inheritParams create_gene_binary
#' @keywords internal
#' @return a data frame
#'
.fusions_gene_binary <- function(fusion,
                                 samples,
                                 specify_panel) {

  fusion <- fusion %>%
    stats::na.omit() %>%
    distinct()

  fus_bm <- .process_binary(data = fusion, samples = samples, type = "fus")

  return(fus_bm)
}


# CNA Binary Matrix -----------------------------------------------------

#' Make Binary Matrix From CNA data frame
#'
#' @inheritParams create_gene_binary
#' @keywords internal
#' @return a data frame
#'
.cna_gene_binary <- function(cna,
                             samples,
                             specify_panel,
                             high_level_cna_only) {

  # * Remove lower level CNA if specified ----
  if (high_level_cna_only) {
    cna2 <- cna %>%
      filter(!(.data$alteration %in% c("loss", "gain") |
        is.na(.data$alteration)))
  } else {
    cna <- cna %>%
      mutate(
        alteration =
          dplyr::case_when(
            .data$alteration == "gain" ~ "amplification",
            .data$alteration == "loss" ~ "deletion",
            TRUE ~ as.character(.data$alteration)
          )
      )
  }


  cna_del <- .process_binary(
    data = cna,
    samples = samples,
    type = "del"
  )

  cna_amp <- .process_binary(
    data = cna,
    samples = samples,
    type = "amp"
  )

  cna_bm <- full_join(cna_del, cna_amp, by = "sample_id") %>%
    mutate(across(-c("sample_id"),
      .fns = function(x) ifelse(is.na(x), 0, x)
    ))

  return(cna_bm)
}


# internal binary matrix creation code for use in .XXX_gene_binary() functions

#' Make a binary matrix from list of samples and genes
#'
#' @inheritParams
#'

