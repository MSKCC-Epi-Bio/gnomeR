
test_that("both sanitize functions run with no errors", {
  expect_error(sanitize_mutation_input(gnomeR::mut), NA)
  expect_error(sanitize_fusion_input(gnomeR::fusion), NA)
})


test_that("both sanitize functions run without warnings",{
  expect_warning(sanitize_mutation_input(gnomeR::mut), NA)
  expect_warning(sanitize_fusion_input(gnomeR::fusion), NA)
})

test_that("test to see what happens if pass sanitize a vector", {
  expect_error(sanitize_mutation_input(gnomeR::mut %>%
                                           select(-Hugo_Symbol)))
  expect_error(sanitize_fusion_input(gnomeR::fusion %>%
                                         select(-Hugo_Symbol)))
})

