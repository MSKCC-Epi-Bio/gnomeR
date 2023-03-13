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
#'   rename_columns()%>%
#'   reformat_fusion()
#'
#' samples <- unique(sv_long1$sampleId)
#'
#'
#' @import dplyr
#' @import stringr
#' @import tidyr



reformat_fusion <- function(fusions) {


  # Checks ----------
  if (!is.data.frame(fusions)) {
    cli::cli_abort("{.code fusion} must be a data.frame")
  }


  .check_required_cols(fusions, c("fusion", "sample_id", "hugo_symbol"))



  # Clean dataset ------------
  fusions_sep <- suppressMessages(fusions %>%
    mutate(
      # remove leading space in fusion var
      fusion = str_trim(.data$fusion),
      # remove endings of names
      fusion2 = case_when(
        endsWith(.data$fusion, " fusion") ~ gsub(".{7}$", "", .data$fusion),
        endsWith(.data$fusion, "-intragenic") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, "-INTRAGENIC") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, "-INTERGENIC") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, "-intergenic") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, " truncation") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, " rearrangement") ~ gsub(".{14}$", "", .data$fusion),
        endsWith(.data$fusion, " fusion - Archer") ~ gsub(".{16}$", "", .data$fusion),
        endsWith(.data$fusion, " duplication") ~ gsub(".{11}$", "", .data$fusion),
        endsWith(.data$fusion, " rearrangement") ~ gsub(".{13}$", "", .data$fusion),
        endsWith(.data$fusion, " truncation") ~ gsub(".{10}$", "", .data$fusion),
        endsWith(.data$fusion, "EZH2(NM_004456) rearrangement exon 5") ~ gsub(".{36}$", "", .data$fusion),
        endsWith(.data$fusion, " PAX5(NM_016734) rearrangement intron 8") ~ gsub(".{38}$", "", .data$fusion),
        TRUE ~ gsub(".{20}$", "", .data$fusion)
      )
    ) %>%
    # keep the original gene-gene info somewhere in the data frame and don't overwrite
    dplyr::rename(event_info = .data$fusion) %>%
    mutate(
      gene_order = .data$fusion2,
      # rename known fusions with two hyphens in the name
      event_info = case_when(
        startsWith(.data$event_info, "NPHP3-AS1-STAG1") ~ "NPHP3_AS1-STAG1",
        startsWith(.data$event_info, "STAG1-NPHP3-AS1") ~ "STAG1-NPHP3_AS1",
        startsWith(.data$event_info, "RICTOR-NIPBL-AS1") ~ "RICTOR-NIPBL_AS1",
        startsWith(.data$event_info, "PMF1-BGLAP-RIT1") ~ "PMF1_BGLAP-RIT1",
        startsWith(.data$event_info, "U2AF1-MCM3AP-AS1") ~ "U2AF1-MCM3AP_AS1",
        startsWith(.data$event_info, "MCM3AP-AS1-U2AF1") ~ "MCM3AP_AS1-U2AF1",
        startsWith(.data$event_info, "MIR365-2-NF1") ~ "MIR365_2-NF1",
        startsWith(.data$event_info, "FAT1-F11-AS1") ~ "FAT1-F11_AS1",
        startsWith(.data$event_info, "TCEB1-DLGAP1-AS5") ~ "TCEB1-DLGAP1_AS5",
        startsWith(.data$event_info, "TAP2-HLA-DOB") ~ "TAP2-HLA_DOB",
        startsWith(.data$event_info, "BIVM-ERCC5-GLUD1") ~ "BIVM_ERCC5-GLUD1",
        startsWith(.data$event_info, "GLUD1-BIVM-ERCC5") ~ "GLUD1-BIVM_ERCC5",
        startsWith(.data$event_info, "CREBBP-AATK-AS1") ~ "CREBBP-AATK_AS1",
        startsWith(.data$event_info, "CDKN2B-AS1-CDKN2A") ~ "CDKN2B_AS1-CDKN2A",
        startsWith(.data$event_info, "ETV6-PRH1-PRR4") ~ "ETV6-PRH1_PRR4",
        startsWith(.data$event_info, "KMT2A-TMPRSS4-AS1") ~ "KMT2A-TMPRSS4_AS1",
        startsWith(.data$event_info, "TMPRSS4-AS1-KMT2A") ~ "TMPRSS4_AS1-KMT2A",
        startsWith(.data$event_info, "RNU6-19P-ETV6") ~ "RNU6_19P-ETV6",
        startsWith(.data$event_info, "RAD21-C8orf37-AS1") ~ "RAD21-C8orf37_AS1",
        startsWith(.data$event_info, "BRWD1-IT2-TMPRSS2") ~ "BRWD1_IT2-TMPRSS2",
        startsWith(.data$event_info, "TMPRSS2-BRWD1-IT2") ~ "TMPRSS2-BRWD1_IT2",
        startsWith(.data$event_info, "CTD-2151A2.1-DROSHA") ~ "CTD_2151A2.1-DROSHA",
        startsWith(.data$event_info, "FOXA1-NKX2-8") ~ "FOXA1-NKX2_8",
        startsWith(.data$event_info, "BRAF-LINC-PINT") ~ "BRAF-LINC_PINT",
        startsWith(.data$event_info, "MIR365-2-NF1") ~ "MIR365_2-NF1",
        startsWith(.data$event_info, "SOX2-OT-SOX2") ~ "SOX2_OT-SOX2",
        TRUE ~ .data$event_info
      )
    ))

  # find any remaining hugo_symbols that have 2 hyphens
  special_case <- fusions_sep %>%
    filter(str_count(.data$fusion2, "-") > 1) %>%
    arrange(.data$sample_id, .data$fusion2) %>%
    select("sample_id", "fusion2") %>%
    unique()

  invest <- special_case$fusion2

  names(invest) <- rep("!", times = length(special_case$sample_id))

  if (nrow(special_case) > 0) {
    cli::cli_abort(c("Some of your hugo_symbols contain more than one '-' and due to duplicate records, the program cannot identify proper gene names. Check which names should be
                     grouped together, and then replace the '-' in the paired group with '_'. Please investigate the following: ", invest))
  }

  # separate the two genes into their own columns
  fusions_sep1 <- suppressWarnings(fusions_sep %>%
    # select(-hugo_symbol)%>%
    tidyr::separate(.data$fusion2, into = c("site1hugo_symbol", "site2hugo_symbol"), "-") %>%
    select(
      "sample_id", "hugo_symbol", "site1hugo_symbol", "site2hugo_symbol", "event_info",
      "gene_order"
    ) %>%
    # filters out any mistakes because this should never be in the first site position
    filter(!(.data$site1hugo_symbol %in% c("repeat", "insufficient"))) %>%
    unique())


  # There are cases where site 1 and 2 were flipped for a sample_id and listed x2
  # ex: TP53-APC vs APC-TP53 for the same sample would create 4 values when it should be 2


  # make each event a row and count events and number of NA events by sample
  fusions_sep2 <- suppressMessages(fusions_sep1 %>%
    select(-"hugo_symbol") %>%
    tidyr::pivot_longer(!.data$sample_id) %>%
    group_by(.data$sample_id) %>%
    mutate(
      count = n(),
      is_na = as.numeric(is.na(.data$value))
    ) %>%
    unique())

  # count unique events by sample
  fusions_sep3 <- suppressMessages(fusions_sep2 %>%
    select("sample_id", "value") %>%
    unique() %>%
    group_by(.data$sample_id) %>%
    mutate(count_unique = n()) %>%
    select("sample_id", "count_unique") %>%
    unique())


  # sort ids based on how many raw vs unique counts per sample


  # merge unique and raw counts together and compare. If diff is not count_na => problem
  fusions_prob_ids <- suppressMessages(fusions_sep2 %>%
    select("sample_id", "is_na") %>%
    group_by(.data$sample_id) %>%
    summarize(
      count_na = sum(.data$is_na),
      count = n()
    ) %>%
    left_join(fusions_sep3) %>%
    # unique events not equal to observed #, and difference isn't the count_na
    filter(.data$count_unique != .data$count &
             .data$count - .data$count_unique != .data$count_na))

  prob_ids <- fusions_prob_ids$sample_id

  # select ids that are not an issue
  fusions_noprob_ids <- suppressMessages(fusions_sep2 %>%
    select("sample_id", "is_na") %>%
    group_by(.data$sample_id) %>%
    summarize(
      count_na = sum(.data$is_na),
      count = n()
    ) %>%
    left_join(fusions_sep3) %>%
    filter(.data$count_unique == .data$count | .data$count -
             .data$count_unique == .data$count_na))

  noprob_ids <- fusions_noprob_ids$sample_id

  # create dataset with noprob ids

  fusions_noprob <- suppressMessages(fusions_sep1 %>%
    filter(.data$sample_id %in% noprob_ids) %>%
    select(-"hugo_symbol") %>%
    arrange(.data$sample_id, .data$site1hugo_symbol) %>%
    unique())



  # retrieve full fusions for prob ids

  fusions_prob <- suppressMessages(fusions_sep1 %>%
    filter(.data$sample_id %in% prob_ids) %>%
    select(-"hugo_symbol") %>%
    arrange(.data$sample_id, .data$site1hugo_symbol) %>%
    unique())


  # split out samples among the prob ids that include intragenic
  probid_okfusions <- fusions_sep1 %>%
    filter(is.na(.data$site2hugo_symbol))

  if (nrow(probid_okfusions) == 0) {
    intragenic <- "no"
  } else {
    intragenic <- "yes"
  }


  # collect fusions that are two gene fusions (aka non-intra or intergenic ex: gene1-gene2)
  non_intra_samps <- switch(intragenic,
    "yes" = {
      setdiff(fusions_sep1, probid_okfusions)
    },
    "no" = {
      fusions_sep1
    }
  )%>%
    unique() %>%
    select(-"hugo_symbol") %>%
    tidyr::pivot_longer(starts_with("site"))%>%
    arrange(.data$sample_id, .data$value, .data$event_info) %>%
    suppressWarnings() %>%
    suppressMessages()


  # now lets deal with the duplicates in the non_intra_samples dataset

  # check that there are multiple events for each fusion
  check <- non_intra_samps %>%
    group_by(.data$sample_id, .data$value) %>%
    summarize(count = n()) %>%
    suppressWarnings() %>%
    suppressMessages()

  # pivot wider
  non_intra_wide <- non_intra_samps %>%
    # create site2 from event info to check that it is correct
    # some will be listed APC-APC or TP53-TP53 based on event info order
    # and we have to drop these matched pairs
    mutate(site2hugo_symbol = sub(".*-", "", .data$gene_order)) %>%
    select("sample_id", "event_info", "gene_order", "value", "site2hugo_symbol") %>%
    # drop pairs
    filter(.data$value != .data$site2hugo_symbol) %>%
    rename(site1hugo_symbol = .data$value) %>%
    unique() %>%
    suppressWarnings() %>%
    suppressMessages()

  # selecting helper functions
  na_list <- list(c(NA, NA))
  not_all_na <- function(x) any(!is.na(x))
  not_all_null <- function(x) any(!is.null(x))

  # loop through the dataset and pivot wider as many times as needed
  # fusionsset columns should look like this
  # ~sample_id, ~g1_1, ~g2_1, ..., ~g[pairnum]_1, ~g2_1, ~g2_2, ..., ~g[pairnum]_2
  # where gx_y defines x = gene order in fusion and y = fusion id
  # (example: TERT-APC would be g1_1 = TERT, gene g2_1 = APC and if the reverse
  # exists we would see g1_2 = APC, g2_2 = TERT)
  one_row_per_samp <- non_intra_wide %>%
    select(-c("event_info", "site1hugo_symbol", "site2hugo_symbol")) %>%
    mutate(genes = .data$gene_order) %>%
    tidyr::separate(.data$genes, into = c("g1", "g2"), "-") %>%
    group_by(.data$sample_id) %>%
    mutate(test = seq_along(.data$gene_order)) %>%
    # pivot wider again to get rid of duplicates
    pivot_wider(id_cols = .data$sample_id, names_from = .data$test,
                values_from = c(.data$g1, .data$g2)) %>%
    ungroup() %>%
    suppressWarnings() %>%
    suppressMessages()


  # find the number of pairs that there are max
  temp1 <- one_row_per_samp %>%
    dplyr::select(-"sample_id")


  # there have to be enough pairs to run through comparisons
  if (length(colnames(temp1)) > 2) {
    pairnum <- length(colnames(temp1)) / 2
    pair_names <- c(paste0("pair", 1:pairnum))



    # now we want to properly pair the genes together in a list and sort to be
    # alphabetical order (ex: TERT-APC should be APC-TERT) and then only select
    # unique events
    for (x in 1:pairnum) {
      pair <- one_row_per_samp %>%
        select("sample_id", ends_with(as.character(x)))

      pair$temp_pair <- list(vector(mode = "list", length = 2))

      for (y in 1:nrow(pair)) {
        gene1 <- as.character(pair[y, 2])
        gene2 <- as.character(pair[y, 3])

        # here if two gene names exist in a pair we want to alphabetize them
        # with sort(), else empty string or single gene name
        pair$temp_pair[y] <- ifelse(is.na(gene1) & is.na(gene2), c(""),
          ifelse(!is.na(gene1) & is.na(gene2), c(gene1, ""),
            list(sort(c(gene1, gene2)))
          )
        )
      }

      # iterate the pair number and name column
      colnames(pair)[colnames(pair) == "temp_pair"] <- paste0("pair", x)

      # join to overalldataset
      one_row_per_samp <- one_row_per_samp %>%
        left_join(pair) %>%
        suppressWarnings() %>%
        suppressMessages()
    }


    # now check to see if any of these events are repeated with genes flipped as
    # described above. If so, we want to remove them

    # create a shelldataset with sample_ids and the first fusion for each
    # sample_id so we can compare to other pairs
    shell <- one_row_per_samp %>%
      select(-c("pair2":last_col()))


    # loop through the number of pairs 2-pairnum
    for (a in 2:pairnum) {
      pair <- one_row_per_samp %>%
        select(c("sample_id", starts_with("pair")))

      comp <- pair_names[1:(a - 1)]

      # set name of tested column (ex: pair3) to test_pair for ease of loop
      colnames(pair)[colnames(pair) == pair_names[a]] <- "test_pair"


      # loop through all pairs before pair[a] ex(pair1 & pair2 if a = 3)
      # and rename them test_comp
      for (c in comp) {
        colnames(pair)[colnames(pair) == c] <- "test_comp"

        # for each sample indataset, set NA if the testing pair is empty
        # or if it is identical to the comparison pair (ex: TERT-APC == TERT-APC)
        # this ensures that the first time the gene occurs is retained
        # and the second is removed (ex: in this situation pair2 = TERT-APC,
        # pair3 = NA)
        for (b in 1:nrow(pair)) {
          pair$test_pair[b] <- ifelse(pair$test_pair[[b]][1] == "", NA,
            ifelse(identical(pair$test_comp[b], pair$test_pair[b]), NA,
              pair$test_pair[b]
            )
          )
        }

        # reset the column name so that test_comp can be used again in next loop
        colnames(pair)[colnames(pair) == "test_comp"] <- c
      }

      # select only the column we tested in this iteration and sample_ids
      pair <- pair %>%
        select("sample_id", "test_pair")

      # reset column name so that test_pair can be used again in next loop
      colnames(pair)[colnames(pair) == "test_pair"] <- pair_names[a]


      # join the column tested to the shell dataset
      shell <- shell %>%
        left_join(pair) %>%
        suppressWarnings() %>%
        suppressMessages()
    }

    # now that we have the pairs, we can clean things up


    one_row_per_samp <- shell %>%
      select(-c(starts_with("g"))) %>%
      unnest(cols = {{ pair_names }}) %>%
      select_if(not_all_na) %>%
      tidyr::pivot_longer(!.data$sample_id) %>%
      group_by(.data$sample_id, .data$name) %>%
      mutate(gene_num = seq_along(.data$value)) %>%
      ungroup() %>%
      pivot_wider(values_from = .data$value, names_from = .data$gene_num) %>%
      # set hugo_symbol names correctly
      rename(site1hugo_symbol = .data$`1`, site2hugo_symbol = .data$`2`) %>%
      suppressWarnings() %>%
      suppressMessages()
  } else {
    one_row_per_samp <- one_row_per_samp %>%
      select_if(not_all_na) %>%
      mutate(pair = "pair1") %>%
      # set hugo_symbol names correctly
      rename(site1hugo_symbol = .data$g1_1, site2hugo_symbol = .data$g2_1)
  }



  # join the new pairings with the rest of thedataset
  # to get back the event info and the gene-order
  new_format <- one_row_per_samp %>%
    left_join(non_intra_wide) %>%
    stats::na.omit() %>%
    suppressWarnings() %>%
    suppressMessages()

  #
  not_working <- one_row_per_samp %>%
    left_join(non_intra_wide) %>%
    filter(is.na(.data$event_info)) %>%
    mutate(site2 = .data$site1hugo_symbol) %>%
    select(-c("site1hugo_symbol", "gene_order", "event_info")) %>%
    rename(site1 = .data$site2hugo_symbol) %>%
    rename(
      site1hugo_symbol = .data$site1,
      site2hugo_symbol = .data$site2
    ) %>%
    left_join(non_intra_wide) %>%
    stats::na.omit() %>%
    suppressWarnings() %>%
    suppressMessages()

  if (length(colnames(temp1)) > 2) {
    # bind all of the fusions together
    to_merge <- new_format %>%
      rbind(not_working) %>%
      select(-.data$name) %>%
      rbind(probid_okfusions %>% select(-"hugo_symbol")) %>%
      rbind(fusions_noprob)
  } else {
    to_merge <- new_format
  }




  # merge all of the datasets together
  # including intragenic, non-repeated events, and above dataset

  fus2 <- to_merge %>%
    select("sample_id", "site1hugo_symbol", "site2hugo_symbol", "event_info") %>%
    unique() %>%
    mutate(across(
      !.data$sample_id,
      ~ str_replace(., "_", "-")
    ))%>%
    rename(fusion = .data$event_info)


  return(fus2)
}
