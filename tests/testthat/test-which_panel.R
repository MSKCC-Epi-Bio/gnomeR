
test_that("test one gene only in IMPACT505", {
  hugo_symbol <- "ZNRF3"
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 5)
  expect_equal(nrow(test), length(hugo_symbol))

})

test_that("test two genes gene only in IMPACT505", {
  hugo_symbol <- c("ZNRF3", "TP53")
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 5)
  expect_equal(nrow(test), length(hugo_symbol))

})
