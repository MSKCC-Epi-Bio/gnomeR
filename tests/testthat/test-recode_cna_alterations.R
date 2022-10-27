cna <- rename_columns(gnomeR::cna)

cna_no_alt_long <- cna %>%
  select(-alteration)

not_allowed_nums <- cna %>%
  mutate(alteration = 3)

mixed_cna_labels <- cna %>%
  select(-sample_id)%>% #changing alteration values so drop patient ID
  mutate(alteration = c(1, rep("neutral", nrow(cna) - 1)))


test_that("test no error thrown with example dataset", {
  expect_error(recode_cna_alterations(cna), NA)

  alt_types <- names(table(recode_cna_alterations(cna)$alteration))
  expect_true("high level amplification" %in% alt_types)
  expect_true("homozygous deletion" %in% alt_types)
  expect_equal(length(alt_types), 2)
})


test_that("test what happens when api (wide) cna given", {
  expect_error(recode_cna_alterations(reformat_cna(cna)))
})


test_that("test that missing alteration throws error", {
  expect_error(recode_cna_alterations(cna_no_alt_long))
})

test_that("test that numbers < -2 or > 2 throw error", {
  expect_error(recode_cna_alterations(not_allowed_nums))
})

test_that("allowed to have numbers and characters in alteration column", {
  expect_error(recode_cna_alterations(mixed_cna_labels), NA)

  mixed_cna_labels <- recode_cna_alterations(mixed_cna_labels)

  count_alterations <- mixed_cna_labels %>%
    group_by(alteration)%>%
    summarize(count = n())%>%
    as.data.frame()

  types_alterations <- names(table(mixed_cna_labels$alteration))


  expect_equal(types_alterations, c("gain", "neutral"))
  expect_equal(count_alterations[1, 2], 1) #only one gain
  expect_equal(count_alterations[2, 2], nrow(mixed_cna_labels) - 1) #all others neutral
})

