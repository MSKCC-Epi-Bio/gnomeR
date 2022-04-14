
samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
binary_matrix_ex <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                               mut_type = "somatic_only", snp_only = FALSE) %>%
  select(TP53, TP53.Del, APC, RB1, RB1.Del)

test_that("add_pathways function works with default input", {
  binmat <- gnomeR::binary_matrix(mutation = gnomeR::mut,
                                  cna = gnomeR::cna,
                                  fusion = gnomeR::fusion)
  expect_error(p <- add_pathways(binary_matrix = binmat), NA)

  expect_equal(setdiff(names(p), names(binmat)), names(gnomeR::pathways))


})


test_that("function can be piped from binary_matrix()", {
  expect_error(gnomeR::binary_matrix(mutation = gnomeR::mut,
                                     cna = gnomeR::cna,
                                     fusion = gnomeR::fusion) %>%
                 add_pathways(), NA)
})

# pathways -------------------------------------------------------------
test_that("pass specific pathways", {

  expect_error(p <- add_pathways(binary_matrix = binary_matrix_ex,
                                      pathways = c("Notch")), NA)
  expect_equal(c(names(binary_matrix_ex), "Notch"), names(p))

  expect_error(p <- add_pathways(binary_matrix = binary_matrix_ex,
                                 pathways = c("Notch", "Myc")), NA)

  expect_warning(p <- add_pathways(binary_matrix = binary_matrix_ex,
                                 pathways = c("Notch", "no")), "*")

  expect_equal(c(names(binary_matrix_ex), "Notch"), names(p))

})

test_that("pass incorrect pathway", {

  expect_message(cust <- add_pathways(binary_matrix = binary_matrix_ex,
                                      custom_pathways = c("TP53", "APC")))

  #check summed only mutations in path
  expect_equal(sum(cust$custom_pathway),
               sum(binary_matrix_ex$TP53, binary_matrix_ex$APC))

  cust2 <- add_pathways(binary_matrix = binary_matrix_ex,
                        custom_pathways = c("TP53", "APC"),
                        count_pathways_by = "gene")

  #check summed only mutations in path and not all
  expect_true(sum(cust2$custom_pathway) > sum(cust$custom_pathway))

})

# custom pathways -------------------------------------------------------------
test_that("expect warning/results when no CNA/Fusions in custom", {

  expect_message(cust <- add_pathways(binary_matrix = binary_matrix_ex,
                              custom_pathways = c("TP53", "APC")))

  #check summed only mutations in path
  expect_equal(sum(cust$custom_pathway),
               sum(binary_matrix_ex$TP53, binary_matrix_ex$APC))

  cust2 <- add_pathways(binary_matrix = binary_matrix_ex,
                                      custom_pathways = c("TP53", "APC"),
                                      count_pathways_by = "gene")

  #check summed only mutations in path and not all
  expect_true(sum(cust2$custom_pathway) > sum(cust$custom_pathway))

})

test_that("vector custom pathway with NULL pathways", {

  expect_message(cust <- add_pathways(binary_matrix = binary_matrix_ex,
                                      pathways = NULL,
                                      custom_pathways = c("TP53", "APC")))

  expect_equal(c(names(binary_matrix_ex), "custom_pathway"),
               names(cust))
})

test_that("list custom pathway with NULL pathways", {

  expect_message(cust <- add_pathways(binary_matrix = binary_matrix_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        "path1" = c("TP53", "APC"),
                                        "path2" = c("RB1.Del"))), "*")

  expect_equal(c(names(binary_matrix_ex), "path1", "path2"),
               names(cust))
})

test_that("list custom pathway with NULL names", {

  expect_message(cust <- add_pathways(binary_matrix = binary_matrix_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        c("TP53", "APC"),
                                        c("RB1.Del"))), "*")

})

test_that("list custom pathway with NULL names", {

  expect_warning(add_pathways(binary_matrix = binary_matrix_ex,
               pathways = NULL,
               count_pathways_by = "gene",
               custom_pathways = list(
                 c("TP53", "APC"),
                 c("RB1.Del"))))

})
# count_pathways_by ----------------------------------------------------------
test_that("works with count_pathways_by gene or alt ", {

  gene <- add_pathways(binary_matrix = binary_matrix_ex,
                       pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "gene")

  expect_message(alt <- add_pathways(binary_matrix = binary_matrix_ex,
                      pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "alteration"))

  expect_equal(alt$custom_pathway, binary_matrix_ex$TP53)
  expect_equal(sum(gene$custom_pathway), sum(binary_matrix_ex$TP53, binary_matrix_ex$TP53.Del))

})

