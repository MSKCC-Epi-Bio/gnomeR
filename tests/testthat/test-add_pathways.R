# Ideas for tests to be added-----------
#
# 1. Test that if custom_pathway with no suffix passed it will error e.g. c(TP53, APC)
# Test the above when list or vector passed to custom_pathway
# example of list
# custom_pathways = list(
#   c("TP53.all", "APC"),
#   c("RB1.Del", "FGFR3.all"))
#
# 2. Test that passing custom pathway with GENE.all would return same as
# custom pathway with GENE.mut, GENE.Del, GENE.Amp, GENE.fus
#
# 3. Test that a binary matrix with columns that have .mut suffix (e.g. column name TP53.mut) will correctly be
# processed same as columns with no suffix  (e.g. column name TP53)
#
# 4. Make sure all works with both a few default `pathways` AND `custom_pathways` passed at the same time

test_that("add_pathways function works with default input", {

  binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10,],
                                  cna = gnomeR::cna,
                                  fusion = gnomeR::sv[1:10,])

  expect_error(p <- add_pathways(gene_binary = binmat), NA)
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


test_that("vector custom pathway with NULL pathways", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = c("TP53.mut", "APC.Del"))

  expect_equal(c(names(gene_binary_ex), "pathway_custom"),
               names(cust))
})

test_that("list custom pathway with NULL pathways", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        "paththingyyy" = c("TP53.mut", "APC.mut"),
                                        "patheroo" = c("RB1.Del")))

  expect_equal(c(names(gene_binary_ex), "pathway_paththingyyy", "pathway_patheroo"),
               names(cust))
})

test_that("list custom pathway with NULL names", {


  gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10, ],
                                               cna = gnomeR::cna,
                                               fusion = gnomeR::sv[1:10,])

  cust <- add_pathways(gene_binary = gene_binary_ex,
                                      pathways = NULL,
                                      custom_pathways = list(
                                        c("TP53.mut", "APC.mut"),
                                        c("RB1.Del")))


  expect_equal(c(names(gene_binary_ex), "pathway_custom_1", "pathway_custom_2"),
               names(cust))

})


# These are tests of old deprecated arguments. They can be recycled or deleted

# # # count_pathways_by ----------------------------------------------------------
# test_that("same results when only mutations passed with count_pathways_by gene or alt ", {
#
#   mut_valid_sample_ids<- unique(gnomeR::mutations$sampleId)[1:10]
#   gene_binary_ex <- create_gene_binary(sample=mut_valid_sample_ids,
#                                        mutation=gnomeR::mutations)
#
#   gene <- add_pathways(gene_binary = gene_binary_ex,
#                        pathways = NULL,
#                        custom_pathways = c("TP53", "APC"),
#                        count_pathways_by = "gene")
#
#   expect_message(alt <- add_pathways(gene_binary = gene_binary_ex,
#                       pathways = NULL,
#                        custom_pathways = c("TP53", "APC"),
#                        count_pathways_by = "alteration"))
#
#   expect_equal(alt$pathway_custom, gene$pathway_custom)
#
# })
#
# test_that("works with count_pathways_by gene or alt ", {
#
#   mut_valid_sample_ids <- c("P-0002375-T01-IM3", "P-0003541-T01-IM5", "P-0005571-T01-IM5")
#   gene_binary_ex <- create_gene_binary(sample=mut_valid_sample_ids,
#                                        mutation=gnomeR::mutations,
#                                        cna = gnomeR::cna)
#
#   gene <- add_pathways(gene_binary = gene_binary_ex,
#                        pathways = NULL,
#                        custom_pathways = c("TP53", "APC"),
#                        count_pathways_by = "gene")
#
#   expect_message(alt <- add_pathways(gene_binary = gene_binary_ex,
#                                      pathways = NULL,
#                                      custom_pathways = c("TP53", "APC"),
#                                      count_pathways_by = "alteration"))
#
#   expect_message(alt2 <- add_pathways(gene_binary = gene_binary_ex,
#                                      pathways = NULL,
#                                      custom_pathways = c("TP53.Del", "APC"),
#                                      count_pathways_by = "alteration"))
#
#   expect_equal(sum(gene$pathway_custom), sum(gene_binary_ex$TP53, gene_binary_ex$TP53.Del))
#   expect_equal(sum(alt$pathway_custom), sum(gene_binary_ex$TP53))
#   expect_equal(sum(alt2$pathway_custom), sum(gene_binary_ex$TP53.Del))
#
#
# })
#



