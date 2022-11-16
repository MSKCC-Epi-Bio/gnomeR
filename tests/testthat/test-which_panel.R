
test_that("test one gene only in IMPACT505", {
  hugo_symbol <- "ZNRF3"
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 5)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_equal(test$IMPACT505, "yes")
  expect_equal(test$IMPACT468, "no")
  expect_equal(test$IMPACT410, "no")
  expect_equal(test$IMPACT341, "no")

})

test_that("test two genes gene only in IMPACT505", {
  hugo_symbol <- c("ZNRF3", "TP53")
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 5)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_equal(test$IMPACT505[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT468[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT410[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT341[test$genes_in_panel == "TP53"], "yes")

})

test_that("test non-impact gene", {
  hugo_symbol <- c("XXX")
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 5)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_equal(test$IMPACT505, "no")
  expect_equal(test$IMPACT468, "no")
  expect_equal(test$IMPACT410, "no")
  expect_equal(test$IMPACT341, "no")

})
