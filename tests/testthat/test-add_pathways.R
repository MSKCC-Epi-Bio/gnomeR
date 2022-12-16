#
# test_that("add_pathways function works with default input", {
#
#   cna_long <- pivot_cna_longer(gnomeR::cna[1:10,])
#   binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mut[1:10,],
#                                   cna = cna_long,
#                                   fusion = gnomeR::fusion[1:10,])
#
#   expect_error(p <- add_pathways(gene_binary = binmat), NA)
# #  expect_equal(setdiff(names(p), names(binmat)), paste0("pathway_", names(gnomeR::pathways)))
#
#
# })
#
#
# test_that("function can be piped from create_gene_binary()", {
#
#
#   cna_long <- pivot_cna_longer(gnomeR::cna[1:10,])
#
#   expect_error(gnomeR::create_gene_binary(mutation = gnomeR::mut[1:10, ],
#                                      cna = cna_long,
#                                      fusion = gnomeR::fusion[1:10,]) %>%
#                  add_pathways(), NA)
# })
#
# # pathways -------------------------------------------------------------
# test_that("pass specific pathways", {
#
#
#   cna_long <- pivot_cna_longer(gnomeR::cna[1:10,])
#   gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mut[1:10, ],
#                                           cna = cna_long,
#                                           fusion = gnomeR::fusion[1:10,])
#
#   expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
#                                       pathways = c("Notch")), NA)
#
#   expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))
#
#   expect_error(p <- add_pathways(gene_binary = gene_binary_ex,
#                                  pathways = c("Notch", "Myc")), NA)
#
#   expect_warning(p <- add_pathways(gene_binary = gene_binary_ex,
#                                  pathways = c("Notch", "no")), "*")
#
#   expect_equal(c(names(gene_binary_ex), "pathway_Notch"), names(p))
#
# })
#
# test_that("pass incorrect pathway", {
#
#   cna_long <- pivot_cna_longer(gnomeR::cna[1:10,])
#   gene_binary_ex <- gnomeR::create_gene_binary(mutation = gnomeR::mut[1:10, ],
#                                                cna = cna_long,
#                                                fusion = gnomeR::fusion[1:10,])
#
#   # create fake ALK column so counting by genes is different than counting by alterations
#   gene_binary_ex$ALK = gene_binary_ex$KRAS
#
#   expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
#                                       custom_pathways = c("TP53", "ALK")))
#
#   #check summed only mutations in path
#   expect_equal(sum(cust$pathway_custom),
#                sum(gene_binary_ex$TP53, gene_binary_ex$ALK))
#
#   cust2 <- add_pathways(gene_binary = gene_binary_ex,
#                         custom_pathways = c("TP53", "ALK"),
#                         count_pathways_by = "gene")
#
#   #check summed only mutations in path and not all
#   expect_true(sum(cust2$pathway_custom) > sum(cust$pathway_custom))
#
#
#   expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
#                               custom_pathways = c("TP53", "ALK")))
#
#   #check summed only mutations in path
#   expect_equal(sum(cust$pathway_custom),
#                sum(gene_binary_ex$TP53, gene_binary_ex$ALK))
#
#   cust2 <- add_pathways(gene_binary = gene_binary_ex,
#                                       custom_pathways = c("TP53", "ALK"),
#                                       count_pathways_by = "gene")
#
#   #check summed only mutations in path and not all
#   expect_true(sum(cust2$pathway_custom) > sum(cust$pathway_custom))
#
# })
#
# test_that("vector custom pathway with NULL pathways", {
#
#   expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
#                                       pathways = NULL,
#                                       custom_pathways = c("TP53", "APC")))
#
#   expect_equal(c(names(gene_binary_ex), "pathway_custom"),
#                names(cust))
# })
#
# test_that("list custom pathway with NULL pathways", {
#
#   expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
#                                       pathways = NULL,
#                                       custom_pathways = list(
#                                         "path1" = c("TP53", "APC"),
#                                         "path2" = c("RB1.Del"))), "*")
#
#   expect_equal(c(names(gene_binary_ex), "pathway_path1", "pathway_path2"),
#                names(cust))
# })
#
# test_that("list custom pathway with NULL names", {
#
#   expect_message(cust <- add_pathways(gene_binary = gene_binary_ex,
#                                       pathways = NULL,
#                                       custom_pathways = list(
#                                         c("TP53", "APC"),
#                                         c("RB1.Del"))), "*")
#
# })
#
# test_that("list custom pathway with NULL names", {
#
#   expect_warning(add_pathways(gene_binary = gene_binary_ex,
#                pathways = NULL,
#                count_pathways_by = "gene",
#                custom_pathways = list(
#                  c("TP53", "APC"),
#                  c("RB1.Del"))))
#
# })
# # count_pathways_by ----------------------------------------------------------
# test_that("works with count_pathways_by gene or alt ", {
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
#   expect_equal(alt$pathway_custom, gene_binary_ex$TP53)
#   expect_equal(sum(gene$pathway_custom),
#                sum(gene_binary_ex$TP53, gene_binary_ex$TP53.Del))
#
# })
#
