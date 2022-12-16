#
#
# #thoughts so far:
# # still need title case for binary matrix function
#
# cbioportalR::set_cbioportal_db("public")
#
#
# test_that("test rename_columns runs with no errors", {
#
#   expect_error(rename_columns(gnomeR::muations), NA)
#
# })
#
# # do we want it to display what it changed? currently does not.
# test_that("test rename_columns does not have messages", {
#
#   expect_message(rename_columns(gnomeR::muations), NA)
#
# })
#
#
# test_that("test colnames are renamed properly", {
#
#   expect_snapshot(
#     waldo::compare(colnames(gnomeR::muations),
#                    colnames(rename_columns(gnomeR::muations)))
#   )
#
# })
#
# #issue is that this is not title case
# test_that("binary matrix runs with renamed columns without error", {
#
#   expect_error(
#    gnomeR::create_gene_binary(samples = gnomeR::muations$Tumor_Sample_Barcode,
#                               mutation = rename_columns(gnomeR::muations)), NA)
#
# })
#
#
