reformat_fusion <- function(data){


  ##########################################################
  ################## clean data#######################
  ##########################################################

  .required_cols(data, "fusion", "sample_id")

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
    mutate(gene_order = fusion2,
           #rename known fusions with
           event_info = case_when(
               startsWith(event_info, "NPHP3-AS1-STAG1") ~ "NPHP3_AS1-STAG1",
               startsWith(event_info, "STAG1-NPHP3-AS1") ~ "STAG1-NPHP3_AS1",
               startsWith(event_info, "RICTOR-NIPBL-AS1") ~ "RICTOR-NIPBL_AS1",
               startsWith(event_info, "PMF1-BGLAP-RIT1") ~ "PMF1_BGLAP-RIT1",
               startsWith(event_info, "U2AF1-MCM3AP-AS1") ~ "U2AF1-MCM3AP_AS1",
               startsWith(event_info, "MCM3AP-AS1-U2AF1") ~ "MCM3AP_AS1-U2AF1",
               startsWith(event_info, "MIR365-2-NF1") ~ "MIR365_2-NF1",
               startsWith(event_info, "FAT1-F11-AS1") ~ "FAT1-F11_AS1",
               startsWith(event_info, "TCEB1-DLGAP1-AS5") ~ "TCEB1-DLGAP1_AS5",
               startsWith(event_info, "TAP2-HLA-DOB") ~ "TAP2-HLA_DOB",
               startsWith(event_info, "BIVM-ERCC5-GLUD1") ~ "BIVM_ERCC5-GLUD1",
               startsWith(event_info, "GLUD1-BIVM-ERCC5") ~ "GLUD1-BIVM_ERCC5",
               startsWith(event_info, "CREBBP-AATK-AS1") ~ "CREBBP-AATK_AS1",
               startsWith(event_info, "CDKN2B-AS1-CDKN2A") ~ "CDKN2B_AS1-CDKN2A",
               startsWith(event_info, "ETV6-PRH1-PRR4") ~ "ETV6-PRH1_PRR4",
               startsWith(event_info, "KMT2A-TMPRSS4-AS1") ~ "KMT2A-TMPRSS4_AS1",
               startsWith(event_info, "TMPRSS4-AS1-KMT2A") ~ "TMPRSS4_AS1-KMT2A",
               startsWith(event_info, "RNU6-19P-ETV6") ~ "RNU6_19P-ETV6",
               startsWith(event_info, "RAD21-C8orf37-AS1") ~ "RAD21-C8orf37_AS1",
               startsWith(event_info, "BRWD1-IT2-TMPRSS2") ~ "BRWD1_IT2-TMPRSS2",
               startsWith(event_info, "TMPRSS2-BRWD1-IT2") ~ "TMPRSS2-BRWD1_IT2",
               startsWith(event_info, "CTD-2151A2.1-DROSHA") ~ "CTD_2151A2.1-DROSHA",
               startsWith(event_info, "FOXA1-NKX2-8") ~ "FOXA1-NKX2_8",
               startsWith(event_info, "BRAF-LINC-PINT") ~ "BRAF-LINC_PINT",
               startsWith(event_info, "MIR365-2-NF1") ~ "MIR365_2-NF1",
               startsWith(event_info, "SOX2-OT-SOX2") ~ "SOX2_OT-SOX2",
               TRUE ~ event_info)))

  special_case <- data_sep %>%
    filter(str_count(fusion2, "-") > 1)%>%
    arrange(sample_id, fusion2)%>%
    select(sample_id, fusion2)%>%
    unique()

  invest <- special_case$sample_id

  names(invest) <- rep("!", times = length(special_case$sample_id))

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

  data_probid_oksamp1 <- data_try %>%
    filter(endsWith(event_info, "INTRAGENIC") | endsWith(event_info, "intragenic") | endsWith(event_info, "INTERGENIC") | endsWith(event_info, "-intragenic - Archer"))%>%
    unique()%>%
    select(-hugo_symbol)%>%
    pivot_wider(id_cols = c(sample_id, gene_order, event_info))%>%
    suppressWarnings()%>%
    suppressMessages()

  # collect samples that are two gene fusions (ex: gene1-gene2)
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

  # see how many events only occur once
  test_sum2 <- non_intra_samps %>%
    group_by(sample_id, value)%>%
    summarize(count = n())%>%
    filter(count == 1)%>%
    select(-count)%>%
    suppressWarnings()%>%
    suppressMessages()


  # now lets deal with the duplicates in the non_intra_samples dataset

  #check that there are multiple events for each fusion
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

  # selecting helper functions
  na_list <- list(c(NA, NA))
  not_all_na <- function(x) any(!is.na(x))
  not_all_null <- function(x) any(!is.null(x))

  #loop through the dataset and pivot wider as many times as needed
  #dataset columns should look like this
  # ~sample_id, ~g1_1, ~g2_1, ..., ~g[pairnum]_1, ~g2_1, ~g2_2, ..., ~g[pairnum]_2
  # where gx_y defines x = gene order in fusion and y = fusion id
  #(example: TERT-APC would be g1_1 = TERT, gene g2_1 = APC and if the reverse
  # exists we would see g1_2 = APC, g2_2 = TERT)
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



 # now we want to properly pair the genes together in a list and sort to be
 # alphabetical order (ex: TERT-APC should be APC-TERT) and then only select
  # unique events
  for(x in 1:pairnum){
    pair <- get_vector %>%
      select(sample_id, ends_with(as.character(x)))

    pair$temp_pair <- list(vector(mode='list', length=2))

    for(y in 1:nrow(pair)){
      gene1 <- as.character(pair[y,2])
      gene2 <- as.character(pair[y,3])

      # here if two gene names exist in a pair we want to alphabetize them
      # with sort(), else empty string or single gene name
      pair$temp_pair[y] <- ifelse(is.na(gene1) & is.na(gene2), c(""),
                                     ifelse(!is.na(gene1) & is.na(gene2), c(gene1, ""),
                                      list(sort(c(gene1, gene2)))))
    }

    # iterate the pair number and name column
    colnames(pair)[colnames(pair) == "temp_pair"] <- paste0("pair", x)

    #join to overall dataset
    get_vector <- get_vector %>%
      left_join(pair)%>%
      suppressWarnings()%>%
      suppressMessages()

  }


  # now check to see if any of these events are repeated with genes flipped as
  # described above. If so, we want to remove them

  # create a shell dataset with sample_ids and the first fusion for each
  # sample_id so we can compare to other pairs
  shell <- get_vector %>%
    select(-c(pair2:last_col()))


  # loop through the number of pairs 2-pairnum
  for(a in 2:pairnum){
    pair <- get_vector %>%
      select(c(sample_id, starts_with("pair")))

    comp <- pair_names[1:(a-1)]

    # set name of tested column (ex: pair3) to test_pair for ease of loop
    colnames(pair)[colnames(pair) == pair_names[a]] <- "test_pair"


    # loop through all pairs before pair[a] ex(pair1 & pair2 if a = 3)
    # and rename them test_comp
      for(c in comp){
      colnames(pair)[colnames(pair) == c] <- "test_comp"

      # for each sample in dataset, set NA if the testing pair is empty
      # or if it is identical to the comparison pair (ex: TERT-APC == TERT-APC)
      # this ensures that the first time the gene occurs is retained
      # and the second is removed (ex: in this situation pair2 = TERT-APC,
      # pair3 = NA)
          for(b in 1:nrow(pair)){
            pair$test_pair[b] <- ifelse(pair$test_pair[[b]][1] == '', NA,
                    ifelse(identical(pair$test_comp[b], pair$test_pair[b]), NA,
                              pair$test_pair[b]))
          }

      # reset the column name so that test_comp can be used again in next loop
      colnames(pair)[colnames(pair) == "test_comp"] <- c

      }

      # select only the column we tested in this iteration and sample_ids
      pair <- pair %>%
        select(sample_id, test_pair)

      # reset column name so that test_pair can be used again in next loop
      colnames(pair)[colnames(pair) == "test_pair"] <- pair_names[a]


      # join the column tested to the shell dataset
      shell <- shell %>%
        left_join(pair)%>%
        suppressWarnings()%>%
        suppressMessages()

    }

  # now that we have the pairs, we can clean things up
  get_vector <- shell %>%
    select(-c(starts_with("g")))%>%
    unnest(cols = {{pair_names}})%>%
    select_if(not_all_na)%>%
    pivot_longer(!sample_id)%>%
    group_by(sample_id, name)%>%
    mutate(gene_num = seq_along(value))%>%
    ungroup()%>%
    pivot_wider(values_from = value, names_from = gene_num)%>%
    # set hugo_symbol names correctly
    rename(site1hugo_symbol = `1`, site2hugo_symbol = `2`)%>%
    suppressWarnings()%>%
    suppressMessages()

  # join the new pairings with the rest of the dataset
  # to get back the event info and the gene-order
  working <- get_vector%>%
    left_join(non_intra_wide) %>%
    na.omit() %>%
    suppressWarnings()%>%
    suppressMessages()

  #
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

  # bind all of the fusions together
  to_merge <- working %>%
    rbind(not_working)%>%
    select(-name)%>%
    suppressWarnings()%>%
    suppressMessages()



  # merge all of the datasets together
  # including intragenic, non-repeated events, and above dataset

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
