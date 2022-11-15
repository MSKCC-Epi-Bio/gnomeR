


test_that("check results for mutations", {

  samp <- unique(gnomeR::mutations$sampleId)
  test_data <- rename_columns(gnomeR::mutations)

  expect_no_error(proc_mut <- .process_binary(data = test_data,
                                  samples = samp,
                                  type = "mut"))

  # check results
  expect_equal(nrow(proc_mut), length(samp))
  expect_equal(ncol(proc_mut)-1, length(unique(gnomeR::mutations$hugoGeneSymbol)))

  x <- gnomeR::mutations %>%
    select(sampleId, hugoGeneSymbol) %>% distinct()

  res_tab <- x$hugoGeneSymbol %>% table() %>% c()
  res_func <- purrr::map_dbl(proc_mut[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})



test_that("check results for fusions", {

  samp <- unique(gnomeR::sv$sampleId)
  test_data <- rename_columns(gnomeR::sv) %>%
    mutate("hugo_symbol" = site_1_hugo_symbol)

  expect_no_error(proc <- .process_binary(data = test_data,
                                              samples = samp,
                                              type = "fus"))

  # check results
  expect_equal(nrow(proc), length(samp))
  expect_equal(ncol(proc)-1, length(unique(test_data$hugo_symbol)))

  x <- test_data %>%
    select(sample_id, hugo_symbol) %>% distinct()

  res_tab <- paste0(x$hugo_symbol, ".fus") %>% table() %>% c()
  res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})


test_that("check results for deletions", {

  samp <- unique(gnomeR::cna$sampleId)
  test_data <- rename_columns(gnomeR::cna) %>%
    mutate(alteration = recode_cna(alteration))

  expect_no_error(proc <- .process_binary(data = test_data,
                                          samples = samp,
                                          type = "del"))
  x <- test_data %>%
    filter(alteration == "deletion") %>%
    select(sample_id, hugo_symbol) %>% distinct()

  # check results
  expect_equal(nrow(proc), length(unique(x$sample_id)))
  expect_equal(ncol(proc)-1, length(unique(x$hugo_symbol)))


  res_tab <- paste0(x$hugo_symbol, ".Del") %>% table() %>% c()
  res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})

test_that("check results for amplifications", {

  samp <- unique(gnomeR::cna$sampleId)
  test_data <- rename_columns(gnomeR::cna) %>%
    mutate(alteration = recode_cna(alteration))

  expect_no_error(proc <- .process_binary(data = test_data,
                                          samples = samp,
                                          type = "amp"))
  x <- test_data %>%
    filter(alteration == "amplification") %>%
    select(sample_id, hugo_symbol) %>% distinct()

  # check results
  expect_equal(nrow(proc), length(unique(x$sample_id)))
  expect_equal(ncol(proc)-1, length(unique(x$hugo_symbol)))


  res_tab <- paste0(x$hugo_symbol, ".Amp") %>% table() %>% c()
  res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})



test_that("sample passed that has events", {

  samp <- c(unique(gnomeR::mutations$sampleId), "ttt")
  test_data <- rename_columns(gnomeR::mutations)

  expect_no_error(proc <- .process_binary(data = test_data,
                                              samples = samp,
                                              type = "mut"))

  x <- test_data %>%
    select(sample_id, hugo_symbol) %>% distinct()

  # check results
  expect_equal(nrow(proc), length(unique(x$sample_id)))
  expect_equal(ncol(proc)-1, length(unique(x$hugo_symbol)))

  res_tab <- x$hugo_symbol %>% table() %>% c()
  res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})


test_that("no events in any samples", {

    samp <- unique(gnomeR::cna$sampleId)

    test_data <- rename_columns(gnomeR::cna) %>%
      mutate(alteration = recode_cna(alteration)) %>%
      filter(alteration == "deletion")

    expect_no_error(proc <- .process_binary(data = test_data,
                                            samples = samp,
                                            type = "amp"))
    x <- test_data %>%
      filter(alteration == "amplification") %>%
      select(sample_id, hugo_symbol) %>% distinct()

    # check results
    expect_equal(nrow(proc), length(unique(x$sample_id)))
    expect_equal(ncol(proc)-1, length(unique(x$hugo_symbol)))


    res_tab <- paste0(x$hugo_symbol, ".Amp") %>% table() %>% c()
    res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
    res_tab <- res_tab[names(res_func)]

    expect_equal(res_tab, res_func)

  })

test_that("only specified list of samples", {
  test_data <- rename_columns(gnomeR::mutations)
  samp <- unique(gnomeR::mutations$sampleId)[1:10]

  expect_no_error(proc <- .process_binary(data = test_data,
                                          samples = samp,
                                          type = "mut"))

  # expect_equal(
  #   create_gene_binary(sample=mut_valid_sample_ids, mutation=gnomeR::mutations) %>%
  #     nrow(),
  #   length(mut_valid_sample_ids))
})
