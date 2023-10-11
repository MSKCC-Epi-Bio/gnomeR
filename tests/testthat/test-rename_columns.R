

# #issue is that this is not title case
test_that("binary matrix runs with renamed columns without error", {

  mut <- gnomeR::mutations %>%
    dplyr::rename("Sample_ID" = sampleId)

  mut2 <- rename_columns(mut)
  expect_false(any(stringr::str_detect(names(mut2), "Sample_ID")))
  expect_true(any(stringr::str_detect(names(mut2), "sample_id")))

  expect_equal(nrow(mut), nrow(mut2))

})

test_that("test that it returns warning message with input data column names", {

  mut2 <- rename_columns(gnomeR::mutations )
  expect_true(!is.null(attr(mut2, "names_dict")))

  cna <- rename_columns(gnomeR::cna)
  expect_true(!is.null(attr(cna, "names_dict")))

  fus <- rename_columns(gnomeR::sv)
  expect_true(!is.null(attr(fus, "names_dict")))


})




test_that("test that it returns warning message with input data column names", {

  mut2 <- capture_messages(create_gene_binary(mutation = gnomeR::mutations))

  expect_true(str_detect(mut2[2], "mutationStatus"))


})
