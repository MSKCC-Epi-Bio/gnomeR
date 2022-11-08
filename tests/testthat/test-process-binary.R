


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

  samp <- unique(gnomeR::cna)
  test_data <- rename_columns(gnomeR::cna) %>%
    .recode_cna_alterations()

  expect_no_error(proc <- .process_binary(data = test_data,
                                          samples = samp,
                                          type = "del"))
  x <- test_data %>%
    filter(alteration == -2) %>%
    select(sample_id, hugo_symbol) %>% distinct()

  # check results
  expect_equal(nrow(proc), length(x$sample_id))
  expect_equal(ncol(proc)-1, length(unique(test_data$hugo_symbol)))



  res_tab <- paste0(x$hugo_symbol, ".fus") %>% table() %>% c()
  res_func <- purrr::map_dbl(proc[, -1], ~sum(.x))
  res_tab <- res_tab[names(res_func)]

  expect_equal(res_tab, res_func)

})
