#' Rename columns from API results to work with gnomeR functions
#'
#' @param df_to_check a data frame to check and recode names as needed
#'
#' @return a renamed data frame
#' @export
#' @examples
#'
#' rename_columns(df_to_check = gnomeR::mutations)
#' rename_columns(df_to_check = gnomeR::sv)
#'
rename_columns <- function(df_to_check) {

  names_df_long <- gnomeR::names_df %>%
    select(contains("_column_name")) %>%
    tidyr::pivot_longer(-"internal_column_name")


  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector
  temp_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(temp_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(temp_dict),
                               temp_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check %>%
      dplyr::rename(!!temp_dict)
  }
}




# Recode Gene Aliases---------------------

#' Recode Hugo Symbol Column
#'
#' @param genomic_df a gene_binary object
#' @param ... Other things passed
#'
#' @return A dataframe with a recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' mut <- rename_columns(gnomeR::mutations)
#'
#' colnames(mut)
#'
#' colnames(recode_alias(genomic_df = mut))
#'
recode_alias <- function(genomic_df, ...) {

  # get table of gene aliases (internal data)
    alias_table <- gnomeR::impact_alias_table %>%
      dplyr::select("hugo_symbol", "alias")

    # recode aliases
    genomic_df$hugo_symbol_old <- genomic_df$hugo_symbol
    genomic_df$hugo_symbol <- purrr::map_chr(genomic_df$hugo_symbol,
                                             ~resolve_alias(gene_to_check = .x,
                                                            alias_table = alias_table))

    message <- genomic_df %>%
      dplyr::filter(.data$hugo_symbol_old != .data$hugo_symbol) %>%
      dplyr::select("hugo_symbol_old", "hugo_symbol") %>%
      dplyr::distinct()


    if(nrow(message) > 0) {
      vec_recode <- purrr::map2_chr(message$hugo_symbol_old,
                                 message$hugo_symbol,
                                 ~paste0(.x, " recoded to ", .y))

      names(vec_recode) <- rep("!", times = length(vec_recode))

      cli::cli_warn(c(
      "To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your dataframe have been recoded (you can prevent this with {.code recode_aliases = FALSE}):",
      vec_recode))

  }

    genomic_df <- genomic_df %>%
      select(-"hugo_symbol_old")

  return(genomic_df)
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
#' resolve_alias("MLL4", alias_table = impact_alias_table)
#'
resolve_alias <- function(gene_to_check, alias_table) {

  if(gene_to_check %in% alias_table$alias) {

    alias_table %>%
      filter(.data$alias == gene_to_check) %>%
      pull("hugo_symbol") %>%
      first() %>%
      as.character()

  } else {
    as.character(gene_to_check)
  }
}


#' IMPACT Panel Annotation of NA's
#'
#' @param gene_binary a processed gene_binary
#'
#' @return  a data frame iwth NAs inserted for genes not tested for given panel versions
#' @export
#'
#'
specify_impact_panels <- function(gene_binary) {

  gene_panels <- gnomeR::gene_panels

  # create data frame of sample IDs
  sample_panel_pair <- rownames(gene_binary) %>%
    as.data.frame() %>%
    stats::setNames("sample_id")

  any_impact <- sum(stringr::str_detect(sample_panel_pair$sample_id,
                                        "-IM|-IH"))

  # check if there are any impact (based on sample ID)
  switch(any_impact == 0,
         cli::cli_abort("There are no IMPACT samples recognized (based on sample_id). If you wish to annotate
                   missingness based on a panel, please pass a data.frame of sample_ids and corresponding panels."))


  # get which IMPACT panel
  sample_panel_pair <- sample_panel_pair %>%
    mutate(panel_id = case_when(
      stringr::str_detect(.data$sample_id, "-IM3") ~ "IMPACT341",
      stringr::str_detect(.data$sample_id, "-IM5") ~ "IMPACT410",
      stringr::str_detect(.data$sample_id, "-IM6") ~ "IMPACT468",
      stringr::str_detect(.data$sample_id, "-IM7") ~ "IMPACT505",
      TRUE ~ "no"
    ))

  # couldn't detect panel
  unk_impact_panel <- sample_panel_pair %>%
    filter(.data$panel_id == "no")

  if(nrow(unk_impact_panel) > 0) {
    cli::cli_alert("Couldn't infer IMPACT panel version from these sample_ids, therefore no NA panel annotation will be done for these: {unk_impact_panel$sample_id}")
  }

  return(sample_panel_pair)
}

#' Annotate Missing Gene Values According to Specific Panels
#'
#' @param sample_panel_pair a data frame of `sample_id`-`panel_id` pairs specifying panels to use for annotation of each sample
#' @param gene_binary a binary matrix of 0/1 indicating alteration yes/no for each sample
#' @return a gene_binary annotated for missingness
#' @export

annotate_any_panel <- function(sample_panel_pair, gene_binary) {

  # if all "no", leave function
  switch(all(sample_panel_pair$panel_id == "no"),
         return(gene_binary))

  sample_panel_pair_nest <- sample_panel_pair %>%
    group_by(.data$panel_id) %>%
    summarise(samples_in_panel = list(.data$sample_id))

  # pull genes for given panels
  panels_needed <- unique(sample_panel_pair_nest$panel_id)

  # has sample IDs and genes for each panel
  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    left_join(gnomeR::gene_panels, by = c("panel_id" = "gene_panel")) %>%
    select(-"entrez_ids_in_panel")

  user_data_genes <- gsub(".fus|.Del|.Amp|.cna", "", colnames(gene_binary))

  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    mutate(na_genes_raw = purrr::map(.data$genes_in_panel,
                                 ~unique(setdiff(user_data_genes, .x)))) %>%
    mutate(na_genes = purrr::map(.data$na_genes_raw,
                                     ~c(
                                       .x,
                                       paste0(.x, ".fus"),
                                       paste0(.x, ".Del"),
                                       paste0(.x, ".Amp"),
                                       paste0(.x, ".cna")
                                     )))


  annotated_data <- purrr::pmap_df(sample_panel_pair_nest,
                                   annotate_specific_panel,
                                   gene_binary = gene_binary)

  return(annotated_data)
}


#' Utility function  to insert NA's According to Panel
#'
#' @param gene_binary a processed binary matrix
#' @param panel_id name of gene panel
#' @param samples_in_panel samples to be annotated for each panel
#' @param na_genes genes to make NA
#' @param ... other args passed
#'
#' @return an annotated data frame
#' @export
#'
#'
annotate_specific_panel <- function(gene_binary,
                           panel_id,
                           samples_in_panel,
                           na_genes, ...) {

  mut_sub <- gene_binary[samples_in_panel, ]
  mut_sub[,stats::na.omit(match(na_genes, colnames(mut_sub)))] <- NA

  return(mut_sub)

}

# create empty data.frame to hold results -----

#' Create empty data.frame to hold results and fill with indicators
#'
#' @param data a dataset with a particular type of genomic event
#' @param samples_final a list of sample IDs inherited from the dataset or specified
#' @param type the types of genomic event in the data set, can be mut/cna/sv
#'
#' @export





.genbin_matrix <- function(data, samples_final,
                           type = c("reformat_cna", "mut", "cna", "sv")) {


  # create vectors to hold import values
  type_alt <- c("deletion", "amplification") #change here if .recode_cna... is accepted
  suffix <- c(".Del", ".Amp", ".fus")
  list_data <- list()
  list_data_new <- list()

  # set unique vectors for each dataset in list

  if (type == "cna") {
    i <- 1 # start counter
  } else {
    i <- length(type_alt) + 1
    if ("site_1_hugo_symbol" %in% colnames(data)) {
      data <- rename(data, "hugo_symbol" = "site_1_hugo_symbol")
    }
    list_data[[1]] <- data
  }


  # if we want cna, need a while loop, else skip by setting counter higher

  while (i <= length(type_alt)) {
    cna_data <- data %>%
      filter(.data$alteration == type_alt[i]) # filter to one type of alt


    list_data[[i]] <- cna_data
    i <- i + 1
  }

  # should happen twice for cna and once for mut and sv
  a <- 1 # set counter for number of datasets in list_data
  for (x in list_data) {

    # assign HS for each dataset, helps with fusion and cna suffixes
    hugo_syms <- unique(x$hugo_symbol)


    if (type == "reformat_cna") {
      data2 <- as.data.frame(matrix(0L,
                                    ncol = length(samples_final) + 1, #+1 for extra col for hugo_symbol names
                                    nrow = length(hugo_syms)
      ))

      colnames(data2) <- c("Hugo_Symbol", samples)
      data2[, 1] <- hugo_syms
      list_data[[1]] <- data2
    } else { #for all other types of data the HS is the co and samp = row
      data2 <- as.data.frame(matrix(0L,
                                    nrow = length(samples_final),
                                    ncol = length(hugo_syms)
      ))
      colnames(data2) <- hugo_syms
      rownames(data2) <- samples_final
    }



    for (y in samples_final) {
      genes <- x$hugo_symbol[x$sample_id %in% y]
      if (length(genes) != 0) {

        # populate matrix depending on type of data
        if (type == "reformat_cna") {
          # fill with any alteration number -2 to 2
          rows_to_fill <- match(unique(as.character(genes)), data2[, 1])
          cols_to_fill <- match(y, colnames(data2))

          data2[rows_to_fill, cols_to_fill] <- x %>%
            filter(.data$sample_id %in% y) %>%
            select("alteration") %>%
            unlist() %>%
            as.numeric()
        } else {
          # fill with 1 if observed else 0
          rows_to_fill <- match(unique(as.character(y)), rownames(data2))
          cols_to_fill <- match(unique(as.character(genes)), colnames(data2))
          data2[rows_to_fill, cols_to_fill] <- 1
        }
      }

      if (length(list_data) > 1) {
        list_data_new[[a]] <- data2 #store dataset for merging
        if(ncol(data2) > 0){
          colnames(list_data_new[[a]]) <- paste0(colnames(list_data_new[[a]]), suffix[a])
        }
      }else{
        list_data_new[[1]] <- data2
      }

    }
    a <- a + 1 # increase counter
  }

  #add fusion suffix
  if (type == "sv"){
    colnames(list_data_new[[1]]) <- paste0(colnames(list_data_new[[1]]), suffix[i])
  }





  if (length(list_data_new) == 1) {
    return(list_data_new[[1]])
  } else {
    i = 2
    list_data_new <- purrr::compact(list_data_new) #drop null datasets from list
    genbin <- list_data_new[[1]]
    while (i <= length(list_data_new)){
      genbin <- list_data_new[[i]] %>% # join all the datasets together
        cbind(genbin)
      i = i + 1
    }

    return(genbin)
  }

}




