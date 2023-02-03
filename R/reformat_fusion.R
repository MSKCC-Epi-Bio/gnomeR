fus_issue_chk <- function(data){


  ##########################################################
  ################## clean data#######################
  ##########################################################

  data_sep <- suppressMessages(data %>%
    mutate(

    #remove leading space in fusion var
      fusion = str_trim(fusion),
    #remove endings to names
    fusion2 = case_when(
      endsWith(fusion, " fusion") ~ gsub('.{7}$', '', fusion),
      endsWith(fusion,"-intragenic") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion,"-INTRAGENIC") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion,"-INTERGENIC") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion,"-intergenic") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion, " truncation") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion, " rearrangement")  ~ gsub('.{14}$', "", fusion),
      endsWith(fusion, " fusion - Archer") ~ gsub('.{16}$', "", fusion),
      endsWith(fusion, " duplication") ~ gsub('.{11}$', "", fusion),
      endsWith(fusion, " rearrangement") ~ gsub('.{13}$', "", fusion),
      endsWith(fusion, " truncation") ~ gsub('.{10}$', "", fusion),
      endsWith(fusion,  "EZH2(NM_004456) rearrangement exon 5") ~ gsub('.{36}$', "", fusion),
      endsWith(fusion, " PAX5(NM_016734) rearrangement intron 8") ~ gsub('.{38}$', "", fusion),

      TRUE ~ gsub('.{20}$', "", fusion)))%>%
    # keep the original gene-gene info somewhere in the data frame and not overwrite
    rename(event_info = fusion) %>%
    mutate(gene_order = fusion2))

  special_case <- data_sep %>%
    filter(str_count(fusion2, "-") > 1)%>%
    arrange(sample_id, fusion2)%>%
    select(sample_id, fusion2)%>%
    unique()

  invest <- purrr::map2_chr(special_case$sample_id,
                            special_case$fusion2,
                            ~paste0(.x, " fusion ", .y))

  names(invest) <- rep("!", times = length(invest))

  if (nrow(special_case) > 0){
    cli::cli_abort(c("Some of your hugo_symbols contain more than one '-' and due to duplicate records, the program cannot identify proper gene names. Check which names should be
                     grouped together, and then replace the '-' in the paired group with '_'. Please investigate the following: ", invest))
  }

  # separate the two genes into their own columns
  data_sep1 <- suppressWarnings(data_sep %>%
    #select(-hugo_symbol)%>%
    separate(fusion2, into = c("site1hugo_symbol", "site2hugo_symbol"),  "-")%>%
    select(sample_id, hugo_symbol, site1hugo_symbol, site2hugo_symbol, event_info,
           gene_order)%>%
    filter(!(site1hugo_symbol %in% c("repeat", "insufficient")))%>%
    unique())

  # get frequency of gene fusion by hugo_symbol
  # There are cases where site 1 and 2 were flipped for a sample_id and listed x2
  # ex: TP53-APC vs APC-TP53 for the same sample would create 4 events when it should be 2
  data_sep_site1 <- suppressMessages(data_sep1 %>%
    select(hugo_symbol, site1hugo_symbol)%>%
    table()%>%
    as.data.frame()%>%
    filter(as.character(hugo_symbol) == as.character(site1hugo_symbol))%>%
    select(hugo_symbol, Freq)%>%
    rename(Freq_site1 = Freq))


  # get frequency for second site and merge to first
  data_sep_site_freq <- suppressMessages(data_sep1 %>%
    select(hugo_symbol, site2hugo_symbol)%>%
    table()%>%
    as.data.frame()%>%
    filter(as.character(hugo_symbol) == as.character(site2hugo_symbol))%>%
    select(hugo_symbol, Freq)%>%
    rename(Freq_site2 = Freq)%>%
    right_join(data_sep_site1))


  # make each event a row and count events and number of NA events by sample
  data_sep2 <- suppressMessages(data_sep1 %>%
    select(-hugo_symbol)%>%
    pivot_longer(!sample_id)%>%
    group_by(sample_id)%>%
    mutate(count = n(),
           is_na = as.numeric(is.na(value)))%>%
    unique())

  # count unique events by sample
  data_sep3 <- suppressMessages(data_sep2 %>%
    select(sample_id, value)%>%
    unique()%>%
    group_by(sample_id)%>%
    mutate(count_unique = n())%>%
    select(sample_id, count_unique)%>%
    unique())





  # sort ids based on how many raw vs unique counts per sample


  # merge unique and raw counts together and compare. If diff is not count_na = prob
  data_prob_ids <- suppressMessages(data_sep2 %>%
    select(sample_id, is_na)%>%
    group_by(sample_id)%>%
    summarize(count_na = sum(is_na),
              count = n())%>%
    left_join(data_sep3)%>%
    filter(count_unique != count & count - count_unique != count_na))

  prob_ids <- data_prob_ids$sample_id

  # select ids that are not an issue
  data_noprob_ids <- suppressMessages(data_sep2 %>%
    select(sample_id, is_na)%>%
    group_by(sample_id)%>%
    summarize(count_na = sum(is_na),
              count = n())%>%
    left_join(data_sep3)%>%
    filter(count_unique == count | count - count_unique == count_na))

  noprob_ids <- data_noprob_ids$sample_id


  # create dataset with noprob ids

  data_noprob <- suppressMessages(data_sep1 %>%
    filter(sample_id %in% noprob_ids)%>%
    select(-hugo_symbol) %>%
    arrange(sample_id, site1hugo_symbol)%>%
    unique())



  # retreive full data for prob ids

  data_prob <- suppressMessages(data_sep1 %>%
    filter(sample_id %in% prob_ids)%>%
    select(-hugo_symbol) %>%
    arrange(sample_id, site1hugo_symbol)%>%
    unique())

  num_dups <- suppressMessages(data_prob %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L))

  data_try <- suppressMessages(data_sep1 %>%
    pivot_longer(starts_with("site"))%>%
    unique()%>%
    arrange(sample_id, gene_order))%>%
    suppressWarnings()


  #split out samples among the prob ids that include intragenic and keep sep

  # 198 events
  data_probid_oksamp1 <- data_try %>%
    filter(endsWith(event_info, "INTRAGENIC") | endsWith(event_info, "intragenic") | endsWith(event_info, "INTERGENIC") | endsWith(event_info, "-intragenic - Archer"))%>%
    unique()%>%
    select(-hugo_symbol)%>%
    pivot_wider(id_cols = c(sample_id, gene_order, event_info))%>%
    suppressWarnings()%>%
    suppressMessages()

  # 404 events, need to check which samples and genes occur x2
  non_intra_samps <- data_try %>%
    filter(!endsWith(event_info, "INTRAGENIC") & !endsWith(event_info, "INTERGENIC"),
           !endsWith(event_info, "intragenic") & !endsWith(event_info, "-intragenic - Archer"))%>%
    unique()%>%
    select(-hugo_symbol)%>%
    arrange(sample_id, value, event_info)%>%
    suppressWarnings()%>%
    suppressMessages()


  # grab non-repeated event with id that has some duplicates
  test_sum <- non_intra_samps %>%
    group_by(sample_id, value)%>%
    summarize(count = n())%>%
    filter(count > 1)%>%
    suppressWarnings()%>%
    suppressMessages()

  test_sum2 <- non_intra_samps %>%
    group_by(sample_id, value)%>%
    summarize(count = n())%>%
    filter(count == 1)%>%
    select(-count)%>%
    suppressWarnings()%>%
    suppressMessages()


  # now lets deal with the duplicates in the non_intra_samples dataset

  #yes, all have multiple events for each gene
  check <- non_intra_samps %>%
    group_by(sample_id, value)%>%
    summarize(count = n())%>%
    suppressWarnings()%>%
    suppressMessages()

  # pivot wider
  non_intra_wide <- non_intra_samps %>%
    mutate(site2hugo_symbol = sub('.*-', '', gene_order))%>%
    select(sample_id, event_info, gene_order, value, site2hugo_symbol)%>%
    filter(value != site2hugo_symbol)%>%
    rename(site1hugo_symbol = value)%>%
    unique()%>%
    suppressWarnings()%>%
    suppressMessages()

  na_list <- list(c(NA, NA))
  not_all_na <- function(x) any(!is.na(x))
  not_all_null <- function(x) any(!is.null(x))



  get_vector <- non_intra_wide %>%
    select(-c(event_info, site1hugo_symbol, site2hugo_symbol))%>%
    mutate(genes = gene_order)%>%
    separate(genes, into = c("g1", "g2"),  "-")%>%
    group_by(sample_id)%>%
    mutate(test = seq_along(gene_order))%>%
    # pivot wider again to get rid of duplicates
    pivot_wider(id_cols = sample_id, names_from = test, values_from = c(g1, g2))%>%
    ungroup()%>%
    suppressWarnings()%>%
    suppressMessages()


  # find the number of pairs that there are max
  temp1 <- get_vector %>%
    dplyr::select(-sample_id)

  pairnum <- length(colnames(temp1))/2
  pair_names <- c(paste0("pair", 1:pairnum))


  for(x in 1:pairnum){
    pair <- get_vector %>%
      select(sample_id, ends_with(as.character(x)))

    pair$temp_pair <- list(vector(mode='list', length=2))

    for(y in 1:nrow(pair)){
      gene1 <- as.character(pair[y,2])
      gene2 <- as.character(pair[y,3])
      pair$temp_pair[y] <- ifelse(is.na(gene1) & is.na(gene2), c(""),
                                     ifelse(!is.na(gene1) & is.na(gene2), c(gene1, ""),
                                      list(sort(c(gene1, gene2)))))
    }

    colnames(pair)[colnames(pair) == "temp_pair"] <- paste0("pair", x)

    get_vector <- get_vector %>%
      left_join(pair)%>%
      suppressWarnings()%>%
      suppressMessages()

  }


  # now check to see if any of these events are repeated with genes flipped
  shell <- get_vector %>%
    select(-c(pair2:last_col()))


  for(a in 2:pairnum){
    pair <- get_vector %>%
      select(c(sample_id, starts_with("pair")))

    comp <- pair_names[1:(a-1)]

    colnames(pair)[colnames(pair) == pair_names[a]] <- "test_pair"



      for(c in comp){
      colnames(pair)[colnames(pair) == c] <- "test_comp"

          for(b in 1:nrow(pair)){
            pair$test_pair[b] <- ifelse(pair$test_pair[[b]][1] == '', NA,
                    ifelse(identical(pair$test_comp[b], pair$test_pair[b]), NA,
                              pair$test_pair[b]))
          }

      colnames(pair)[colnames(pair) == "test_comp"] <- c

      }

      pair <- pair %>%
        select(sample_id, test_pair)

      colnames(pair)[colnames(pair) == "test_pair"] <- pair_names[a]



      shell <- shell %>%
        left_join(pair)%>%
        suppressWarnings()%>%
        suppressMessages()

    }

  get_vector <- shell %>%
    select(-c(starts_with("g")))%>%
    unnest(cols = {{pair_names}})%>%
    select_if(not_all_na)%>%
    pivot_longer(!sample_id)%>%
    group_by(sample_id, name)%>%
    mutate(gene_num = seq_along(value))%>%
    ungroup()%>%
    pivot_wider(values_from = value, names_from = gene_num)%>%
    rename(site1hugo_symbol = `1`, site2hugo_symbol = `2`)%>%
    suppressWarnings()%>%
    suppressMessages()


  working <- get_vector%>%
    left_join(non_intra_wide) %>%
    na.omit() %>%
    suppressWarnings()%>%
    suppressMessages()#get back the event info and the gene-order

  not_working <- get_vector%>%
    left_join(non_intra_wide) %>%
    filter(is.na(event_info))%>%
    mutate(site2 = site1hugo_symbol)%>%
    select(-c(site1hugo_symbol, gene_order, event_info))%>%
    rename(site1 = site2hugo_symbol)%>%
    rename(site1hugo_symbol = site1,
           site2hugo_symbol = site2)%>%
    left_join(non_intra_wide)%>%
    na.omit()%>%
    suppressWarnings()%>%
    suppressMessages()

  to_merge <- working %>%
    rbind(not_working)%>%
    select(-name)%>%
    suppressWarnings()%>%
    suppressMessages()



  # merge all of the datasets together

  data_fus2 <- to_merge %>%
    rbind(data_probid_oksamp1)%>%
    rbind(data_noprob)%>%
    arrange(sample_id, site1hugo_symbol, site2hugo_symbol)%>%
    select(-gene_order)%>%
    unique()%>%
    suppressWarnings()%>%
    suppressMessages()


  return(data_fus2)

  }
