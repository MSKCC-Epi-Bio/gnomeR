
samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
gene_binary_ex <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
                               mut_type = "somatic_only", snp_only = FALSE) %>%
  select(TP53, TP53.Del, APC, RB1, RB1.Del)

test_that("add_pathways function works with default input", {
  binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mut,
                                  cna = gnomeR::cna,
                                  fusion = gnomeR::fusion)
  expect_error(p <- add_pathways(gene_binary = binmat), NA)

  expect_equal(setdiff(names(p), names(binmat)), paste0("pathway_", names(gnomeR::pathways)))


})


test_that("function can be piped from create_gene_binary()", {
  expect_error(gnomeR::create_gene_binary(mutation = gnomeR::mut,
                                     cna = gnomeR::cna,
                                     fusion = gnomeR::fusion) %>%
                 add_pathways(), NA)
})

# pathways -------------------------------------------------------------
test_that("pass specific pathways", {

  expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = c("Notch")), NA)

  expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))

  expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
                                 pathways = c("Notch", "Myc")), NA)

  expect_warning(p <- add_pathways(gene_binary = gene_binary_ex,
                                 pathways = c("Notch", "no")), "*")

  expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))

})

test_that("pass incorrect pathway", {

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      custom_pathways = c("TP53", "APC")))

  #check summed only mutations in path
  expect_equal(sum(cust$pathway_custom),
               sum(gene_binary_ex$TP53, gene_binary_ex$APC))

  cust2 <- add_pathways(gene_binary = gene_binary_ex,
                        custom_pathways = c("TP53", "APC"),
                        count_pathways_by = "gene")

  #check summed only mutations in path and not all
  expect_true(sum(cust2$pathway_custom) > sum(cust$pathway_custom))

})

# custom pathways -------------------------------------------------------------
test_that("expect warning/results when no CNA/Fusions in custom", {

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                              custom_pathways = c("TP53", "APC")))

  #check summed only mutations in path
  expect_equal(sum(cust$pathway_custom),
               sum(gene_binary_ex$TP53, gene_binary_ex$APC))

  cust2 <- add_pathways(gene_binary = gene_binary_ex,
                                      custom_pathways = c("TP53", "APC"),
                                      count_pathways_by = "gene")

  #check summed only mutations in path and not all
  expect_true(sum(cust2$pathway_custom) > sum(cust$pathway_custom))

})

test_that("vector custom pathway with NULL pathways", {

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = c("TP53", "APC")))

  expect_equal(c(names(gene_binary_ex), "pathway_custom"),
               names(cust))
})

test_that("list custom pathway with NULL pathways", {

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        "path1" = c("TP53", "APC"),
                                        "path2" = c("RB1.Del"))), "*")

  expect_equal(c(names(gene_binary_ex), "pathway_path1", "pathway_path2"),
               names(cust))
})

test_that("list custom pathway with NULL names", {

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        c("TP53", "APC"),
                                        c("RB1.Del"))), "*")

})

test_that("list custom pathway with NULL names", {

  expect_warning(add_pathways(gene_binary = gene_binary_ex,
               pathways = NULL,
               count_pathways_by = "gene",
               custom_pathways = list(
                 c("TP53", "APC"),
                 c("RB1.Del"))))

})
# count_pathways_by ----------------------------------------------------------
test_that("works with count_pathways_by gene or alt ", {

  gene <- add_pathways(gene_binary = gene_binary_ex,
                       pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "gene")

  expect_message(alt <- add_pathways(gene_binary = gene_binary_ex,
                      pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "alteration"))

  expect_equal(alt$pathway_custom, gene_binary_ex$TP53)
  expect_equal(sum(gene$pathway_custom),
               sum(gene_binary_ex$TP53, gene_binary_ex$TP53.Del))

})

