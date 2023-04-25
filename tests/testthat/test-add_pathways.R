
test_that("only accecpts tbl_gene_binary objects", {
  fake <- data.frame(sample_id = c(rep("samp", 5)),
                     TERT = c(rep(1, 3), 0, NA))

  expect_error(add_pathways(fake))

  binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10,],
                                       cna = gnomeR::cna,
                                       fusion = gnomeR::sv[1:10,])

  expect_no_error(add_pathways(gene_binary = binmat))

  #test with passing an object in parts

  binmat <- binmat %>%
    summarize_by_gene()

  binmat <- binmat %>%
    subset_by_frequency()

  expect_no_error(test <- add_pathways(gene_binary = binmat))

  expect_true(inherits(test, "tbl_gene_binary"))

})

test_that("produces tbl_gene_binary object", {
  fake <- data.frame(sample_id = c(rep("samp", 5)),
                     TERT = c(rep(1, 3), 0, NA))

  expect_error(add_pathways(fake))

})


test_that("add_pathways function works with default input", {

  # discrepancy between . and - in names in pathway
  binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:20,],
                                  cna = gnomeR::cna,
                                  fusion = gnomeR::sv[1:10,])

  expect_no_error(p <- add_pathways(gene_binary = binmat, count_pathways_by = "alteration"))

  expect_equal(setdiff(names(p), names(binmat)), paste0("pathway_", names(gnomeR::pathways)))


})


test_that("function can be piped from create_gene_binary()", {

  expect_error(gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                     cna = gnomeR::cna,
                                     fusion = gnomeR::sv[1:10,]) %>%
                 add_pathways(), NA)
})

# # pathways -------------------------------------------------------------
test_that("pass specific pathways", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                          cna = gnomeR::cna,
                                          fusion = gnomeR::sv[1:10,])

  expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = c("Notch")), NA)

  expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))

  expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
                                 pathways = c("Notch", "Myc")), NA)

  expect_warning(p <- add_pathways(gene_binary = gene_binary_ex,
                                 pathways = c("Notch", "no")), "Ignoring*")

  expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))

})


test_that("pass incorrect pathway", {

  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:50, ],
                                               cna = cna,
                                               fusion = gnomeR::sv[1:10,])

  # create fake ALK column so counting by genes is different than counting by alterations
  gene_binary_ex$PTEN = gene_binary_ex$AR

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      custom_pathways = c("AKT1", "PTEN")))

  # check summed only mutations in path
  expect_equal(sum(cust$pathway_custom),
               sum(gene_binary_ex$AKT1, gene_binary_ex$PTEN))

  cust2 <- add_pathways(gene_binary = gene_binary_ex,
                        custom_pathways = c("AKT1", "PTEN"),
                        count_pathways_by = "gene")

  #check summed only mutations in path and not all
  expect_true(sum(cust2$pathway_custom) > sum(cust$pathway_custom))


})

test_that("vector custom pathway with NULL pathways", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = c("TP53", "APC")))

  expect_equal(c(names(gene_binary_ex), "pathway_custom"),
               names(cust))
})

test_that("list custom pathway with NULL pathways", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        "paththingyyy" = c("TP53", "APC"),
                                        "patheroo" = c("RB1.Del"))), "*")

  expect_equal(c(names(gene_binary_ex), "pathway_paththingyyy", "pathway_patheroo"),
               names(cust))
})

test_that("list custom pathway with NULL names", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        c("TP53", "APC"),
                                        c("RB1.Del"))), "*")


  expect_equal(c(names(gene_binary_ex), "pathway_custom_1", "pathway_custom_2"),
               names(cust))

})

test_that("list custom pathway with NULL names by gene", {

  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  expect_warning(add_pathways(gene_binary = gene_binary_ex,
               pathways = NULL,
               count_pathways_by = "gene",
               custom_pathways = list(
                 c("TP53", "APC"),
                 c("RB1.Del"))))

})

# # count_pathways_by ----------------------------------------------------------
test_that("same results when only mutations passed with count_pathways_by gene or alt ", {

  mut_valid_sample_ids<- unique(gnomeR::mutations$sampleId)[1:10]
  gene_binary_ex <- create_gene_binary(sample=mut_valid_sample_ids,
                                       mutation=gnomeR::mutations)

  gene <- add_pathways(gene_binary = gene_binary_ex,
                       pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "gene")

  expect_message(alt <- add_pathways(gene_binary = gene_binary_ex,
                      pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "alteration"))

  expect_equal(alt$pathway_custom, gene$pathway_custom)

})

test_that("works with count_pathways_by gene or alt ", {

  mut_valid_sample_ids <- c("P-0002375-T01-IM3", "P-0003541-T01-IM5", "P-0005571-T01-IM5")
  gene_binary_ex <- create_gene_binary(sample=mut_valid_sample_ids,
                                       mutation=gnomeR::mutations,
                                       cna = gnomeR::cna)

  gene <- add_pathways(gene_binary = gene_binary_ex,
                       pathways = NULL,
                       custom_pathways = c("TP53", "APC"),
                       count_pathways_by = "gene")

  expect_message(alt <- add_pathways(gene_binary = gene_binary_ex,
                                     pathways = NULL,
                                     custom_pathways = c("TP53", "APC"),
                                     count_pathways_by = "alteration"))

  expect_message(alt2 <- add_pathways(gene_binary = gene_binary_ex,
                                     pathways = NULL,
                                     custom_pathways = c("TP53.Del", "APC"),
                                     count_pathways_by = "alteration"))

  expect_equal(sum(gene$pathway_custom), sum(gene_binary_ex$TP53, gene_binary_ex$TP53.Del))
  expect_equal(sum(alt$pathway_custom), sum(gene_binary_ex$TP53))
  expect_equal(sum(alt2$pathway_custom), sum(gene_binary_ex$TP53.Del))


})



