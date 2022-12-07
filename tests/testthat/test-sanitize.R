
test_that("both sanitize functions run with no errors", {
  expect_error(sanitize_mutation_input(gnomeR::mutations, include_silent = FALSE), NA)
  expect_error(sanitize_mutation_input(gnomeR::mutations, include_silent = TRUE), NA)
  expect_error(sanitize_fusion_input(gnomeR::sv), NA)
})


test_that("both sanitize functions run without warnings",{
  expect_warning(sanitize_mutation_input(gnomeR::mutations, include_silent = FALSE), NA)
  expect_warning(sanitize_mutation_input(gnomeR::mutations, include_silent = TRUE), NA)
  expect_warning(sanitize_fusion_input(gnomeR::sv), NA)
})

test_that("test to see what happens if pass sanitize a vector", {
  expect_error(sanitize_mutation_input(gnomeR::mutations %>%
                                           select(-Hugo_Symbol),
                                       include_silent = FALSE))
  expect_error(sanitize_fusion_input(gnomeR::sv %>%
                                         select(-Hugo_Symbol)))
})

test_that("alterations properly recoded using internal func", {
  cna <- gnomeR::cna[1:10,]%>%
    sanitize_cna_input()

  expect_true("amplification" %in% names(table(cna$alteration)))
  expect_true("deletion" %in% names(table(cna$alteration)))

  table <- as.data.frame(table(cna$alteration))

  expect_equal(table$Freq[table$Var1 == "deletion"], 7)
  expect_equal(table$Freq[table$Var1 == "amplification"], 3)
})


