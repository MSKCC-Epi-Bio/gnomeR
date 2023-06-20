
test_that("test one gene only in IMPACT505", {
  hugo_symbol <- "ZNRF3"
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 7)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_equal(test$IMPACT505, "yes")
  expect_equal(test$IMPACT468, "no")
  expect_equal(test$IMPACT410, "no")
  expect_equal(test$IMPACT341, "no")
  expect_equal(test$`IMPACT-HEME-400`, "no")
  expect_equal(test$`IMPACT-HEME-468`, "no")
})

test_that("test vector of genes including recodes", {
  hugo_symbol <- c("SOS1", "AGO2", "MLL3", "MLL", "TP53")

  expect_message(test <- which_impact_panel(hugo_symbol))

  expect_equal(ncol(test), 7)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_false("MLL3" %in% test$genes_in_panel)
  expect_false("MLL" %in% test$genes_in_panel)
  expect_true("KMT2C" %in% test$genes_in_panel)

  expect_equal(test$IMPACT505[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT468[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT410[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$IMPACT341[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$`IMPACT-HEME-400`[test$genes_in_panel == "TP53"], "yes")
  expect_equal(test$`IMPACT-HEME-468`[test$genes_in_panel == "TP53"], "yes")

})

test_that("test non-impact gene", {
  hugo_symbol <- c("XXX")
  test <- which_impact_panel(hugo_symbol)

  expect_equal(ncol(test), 7)
  expect_equal(nrow(test), length(hugo_symbol))

  expect_equal(test$IMPACT505, "no")
  expect_equal(test$IMPACT468, "no")
  expect_equal(test$IMPACT410, "no")
  expect_equal(test$IMPACT341, "no")
  expect_equal(test$`IMPACT-HEME-400`, "no")
  expect_equal(test$`IMPACT-HEME-468`, "no")
})
