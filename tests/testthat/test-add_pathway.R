test_that("add_pathways function works with basic input", {
  binmat <- gnomeR::binary_matrix(mutation = gnomeR::mut,
                                  cna = gnomeR::cna,
                                  fusion = gnomeR::fusion)
  expect_error(add_pathway(bin_mat = binmat, pathway = "all"), NA)
})


test_that("function can be piped from binary_matrix()", {
  expect_error(gnomeR::binary_matrix(mutation = gnomeR::mut,
                                     cna = gnomeR::cna,
                                     fusion = gnomeR::fusion) %>%
                 add_pathway(pathway = "all"), NA)
})


test_that("warning is printed when all is included with other pathways", {
  expect_warning(add_pathway(bin_mat = binmat, pathway = c("RTK/RAS", "Notch", "all")))
})


test_that("function works with only mutation data", {
  expect_error(gnomeR::binary_matrix(mutation = gnomeR::mut) %>%
                 add_pathway(bin_mat = binmat, pathway = c("all")))
})
