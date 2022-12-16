



test_that("test no error thrown with example dataset", {
  cna <- rename_columns(gnomeR::cna[1:10, ]) %>%
    mutate(alteration = as.character(alteration))

  recoded_cna <- expect_error(mutate(cna, alteration = recode_cna(alteration)), NA)

  alt_types <- names(table(recoded_cna$alteration))
  expect_true("amplification" %in% alt_types)
  expect_true("deletion" %in% alt_types)
  expect_equal(length(alt_types), 2)

})


test_that("test that numbers < -2 or > 2 throw error", {

  not_allowed_nums <- gnomeR::cna[1:10, ] %>%
    mutate(alteration = "3")

  expect_error(mutate(not_allowed_nums, alteration = recode_cna(alteration)), "Unknown*")
})

test_that("allowed to have numbers and characters in alteration column", {

  mixed_cna_labels <- gnomeR::cna[1:10, ] %>%
    rename_columns(.) %>%
    select(-sample_id) %>% #changing alteration values so drop patient ID
    mutate(alteration = c(1, -2, rep("neutral", nrow(.) - 2)))


  t1 <- table(mixed_cna_labels$alteration)

  expect_error(mixed_relabeled <- mutate(mixed_cna_labels, alteration = recode_cna(alteration)), NA)

  t2 <- table(mixed_relabeled$alteration)

  expect_true(length(setdiff(t1, t2)) == 0)

  types_alterations <- names(table(mixed_cna_labels$alteration))


})


test_that("missing data in CNA is retained in results", {

  na_cna <- gnomeR::cna[1:10, ] %>%
    rename_columns(.) %>%
    select(-sample_id) %>% #changing alteration values so drop patient ID
    mutate(alteration = c(1, NA, rep("neutral", nrow(.) - 2)))


  expect_no_error(nas_data <- mutate(na_cna, alteration = recode_cna(alteration)))
  expect_equal(sum(is.na(nas_data$alteration)), 1)

  unk_cna <- gnomeR::cna[1:10, ] %>%
    rename_columns(.) %>%
    select(-sample_id) %>% #changing alteration values so drop patient ID
    mutate(alteration = c(1, "unknown", rep("neutral", nrow(.) - 2)))

  expect_error(mutate(unk_cna, alteration = recode_cna(alteration)))

})
