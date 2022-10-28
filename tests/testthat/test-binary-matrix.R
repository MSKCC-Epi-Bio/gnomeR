#
# # Test Binary Matrix Arguments -----------------------------------------------------------
#
# # General tests ---
# test_that("test gene_binary with mutation runs with no errors", {
#
#   expect_error(create_gene_binary(mutation = gnomeR::cbp_mut, recode_aliases = FALSE), NA)
#
#   res_mut <- create_gene_binary(mutation = gnomeR::cbp_mut, recode_aliases = FALSE)
#   expect_true(length(res_mut) > 0)
# })
#
#
# test_that("test gene_binary with cna runs with no errors", {
#
#
#   expect_error(create_gene_binary(cna = gnomeR::cbp_cna, recode_aliases = FALSE), NA)
#
#   res <- create_gene_binary(cna = gnomeR::cbp_cna, recode_aliases = FALSE)
#   expect_true(length(res) > 0)
# })
#
#
# test_that("test gene_binary with fusions runs with no errors", {
#
#   expect_error(create_gene_binary(fusion = gnomeR::cbp_sv, recode_aliases = FALSE), NA)
#
#   res_fusion <- create_gene_binary(fusion = gnomeR::fusion, recode_aliases = FALSE)
#   expect_true(length(res_fusion) > 0)
#
# })
#
# test_that("check cna with no alterations are omitted from results", {
#
#   res <- create_gene_binary(mutation = gnomeR::mut,
#                        cna = gnomeR::cna,
#                        fusion = gnomeR::fusion)
#   cna_ids <- names(gnomeR::cna)[-1] %>%
#     str_replace_all(fixed("."), "-")
#
#   omitted_ids <- setdiff(cna_ids, rownames(res))
#
#   omitted_ids <- omitted_ids %>%
#     str_replace_all(fixed("-"), fixed("."))
#
#   check_they_are_zero <- gnomeR::cna %>% select(all_of(omitted_ids)) %>%
#     purrr::map_dbl(., ~sum(.x))
#
#   expect_true(sum(check_they_are_zero) == 0)
# })
#
# # test samples argument ----
# # what happens when you pass a vector? What about if you don't specify it (don't pass anything)?
# # what happens when you pass impact samples (-IM5/IM6/IM7 etc)?  non impact samples? A mix?
#
#
# test_that("Check create_gene_binary() provide specific sample data if pass a vector", {
#
#   #what happens when you pass a vector?
#   mut_valid_sample_ids<-create_gene_binary( mutation= gnomeR::mut) %>%
#                      rownames() %>%
#                      head(n=10)
#
#   expect_equal(
#     create_gene_binary(sample=mut_valid_sample_ids, mutation=gnomeR::mut) %>%
#       nrow(),
#     length(mut_valid_sample_ids))
#
#   #what about if you don't specify it (don't pass anything)?
#
#   expect_lte(
#     create_gene_binary(mutation=gnomeR::mut) %>%
#       nrow(),
#     length(gnomeR::mut[['Tumor_Sample_Barcode']]))
#
# })
#
# test_that("Check create_gene_binary() if sample entered in `sampl_id` with zero mutations/fusions/cna", {
#
#   #what happens when you pass a vector?
#   mut_valid_sample_ids<-create_gene_binary( mutation= gnomeR::mut) %>%
#     rownames() %>%
#     head(n=10)
#
#   add_no_mut_sample <- c(mut_valid_sample_ids, "no_mutations_fake_sample")
#   gene_binary_with_zero <-  create_gene_binary(sample=add_no_mut_sample, mutation=gnomeR::mut)
#
#   sum(gene_binary_with_zero[nrow(gene_binary_with_zero), ])
#   expect_equal(
#     sum(gene_binary_with_zero[nrow(gene_binary_with_zero), ]), 0)
#
#   # should be one more obs in data frame with samples arg specified
#   expect_equal(nrow(gene_binary_with_zero) -1, length(mut_valid_sample_ids))
#
#   # with no fusions in select sample---------
#   gene_binary_with_zero <-  create_gene_binary(samples=add_no_mut_sample,
#                                             mutation = gnomeR::mut,
#                                             cna = gnomeR::cna,
#                                             fusion = gnomeR::fusion)
#
#   expect_false(any(str_detect(names(gene_binary_with_zero), ".fus")))
#   expect_equal(nrow(gene_binary_with_zero), length(add_no_mut_sample))
#
#
#   # with no cna in select sample---------
#   cna_samp <- cna[, c(1, 100)]
#   gene_binary_with_zero <-  create_gene_binary(samples=add_no_mut_sample,
#                                             mutation = gnomeR::mut,
#                                             cna = cna_samp,
#                                             fusion = gnomeR::fusion)
#   expect_false(any(str_detect(names(gene_binary_with_zero), ".Amp")))
#   expect_false(any(str_detect(names(gene_binary_with_zero), "Del")))
#   expect_false(any(str_detect(names(gene_binary_with_zero), ".cna")))
#   expect_equal(nrow(gene_binary_with_zero), length(add_no_mut_sample))
#
# })
#
# # test with and without mut/fusion/cna args passed ----
# # Functions should work with any one of the three passed
# # Does it return results as expected?
# # Trying with and without passing samples arg as well- does it return what you'd expect?
# # what if mut is passed but doesn't have any rows? no columns? no rows or cols?
# test_that("test inputting mut/fusion/cna args can leads to a data.frame output", {
#
#   #Can we obtaine correct result format when either mut/fusion/cna inputted
#   expect_error(create_gene_binary())
#
#   expect_true( create_gene_binary( mutation = gnomeR::mut) %>% is.data.frame())
#
#   expect_true( create_gene_binary( fusion = gnomeR::fusion) %>% is.data.frame())
#
#   expect_true( create_gene_binary( cna = gnomeR::cna) %>% is.data.frame())
#
#   # What if there is no row/col in the file passing to mutation
#
#   #note: there is no error message if a inputting mut data is 0 rows;
#   #      it only return a 0 row and 0 column result
#   expect_equal( gnomeR::mut %>%
#                   dplyr::filter(.data$Hugo_Symbol=="XXXXXXXXX") %>%
#                   create_gene_binary(mutation = .) %>%
#                   nrow(), 0)
#
#   expect_error( gnomeR::mut %>%
#                 dplyr::select(NULL) %>%
#                 create_gene_binary(mutation = .))
#
#
# })
#
# # test mut_type argument ----
#
# test_that("test incorrectly specified arg", {
#
#   expect_error(create_gene_binary(mutation = mut2,mut_type = "somatic_only",
#                              specify_panel = "no"))
# })
#
#
# test_that("test inclusion of NAs in mut_type ", {
#   mut2 = gnomeR::mut
#   mut2$Mutation_Status[1:10]<-NA
#   mut2$Mutation_Status[11:15]<-""
#
#   #example test
#   expect_warning(create_gene_binary(mutation = mut2, specify_panel = "no"), "15 mutations*")
# })
#
#
#
# test_that("test inclusion of NAs in mut_type ", {
#
#   mut2 = gnomeR::mut
#   mut2$Mutation_Status[1:10]<-NA
#   mut2$Mutation_Status[11:15]<-""
#
#   # NA included by default (germline_omitted)
#   expect_warning(see <- create_gene_binary(mutation = mut2, specify_panel = "no"))
#   check <-see$TP53[which(rownames(see)=="P-0000062-T01-IM3")]
#   expect_equal(check, 1)
#
# })
#
# test_that("test inclusion of NAs in mut_type ", {
#   mut2 = gnomeR::mut
#   mut2$Mutation_Status[1:10]<-NA
#   mut2$Mutation_Status[11:15]<-""
#
#
#   # NA included with all
#   see = create_gene_binary(mutation = mut2, specify_panel = "no", mut_type = "all")
#   expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],1)
#
#
#   # NA no longer included with somatic_only
#   see = create_gene_binary(mutation = mut2, mut_type = "somatic_only", specify_panel = "no")
#   expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],0)
#
#   # NA no longer included with germline_only
#   see = create_gene_binary(mutation = mut2, mut_type = "germline_only", specify_panel = "no")
#   expect_equal(ncol(see), 0)
#
#
# })
#
# # test snp_only arg----
# # add general tests
# # What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs
# test_that("test the snp_only arg", {
#
#   #general tests: input T or F (default is F)
#   expect_error( create_gene_binary(mutation=gnomeR::mut, snp_only = T), NA)
#
#   expect_warning( create_gene_binary(mutation=gnomeR::mut,
#                                      snp_only = T,
#                                      recode_aliases = FALSE), NA)
#
#   #0 col return if 0 SNP been inputted or there is no SNP category in Variat Type been inputted
#   mut_snp_zero<- gnomeR::mut %>%
#                  dplyr::filter(.data$Variant_Type!="SNP")
#
#   mut_snp_na<-mut_snp_zero
#   mut_snp_na$Variant_Type<- droplevels(mut_snp_na$Variant_Type)
#
#   expect_equal( create_gene_binary(mutation=mut_snp_zero, snp_only = T) %>%
#                 ncol(),
#                 0)
#
#   expect_equal( create_gene_binary(mutation=mut_snp_na, snp_only = T) %>%
#                   ncol(),
#                 0)
#
#   #What if NA for Variant Type?
#   # note: without Variant Type, the create_gene_binary() still run without error
#   #       snp_only=F will provide full list results and snp_only=T will provide 0 col result
#
#   mut_vt_na<- gnomeR::mut %>%
#               dplyr::mutate(Variant_Type=NA)
#
#   expect_true( create_gene_binary(mutation = mut_vt_na, snp_only = F) %>%
#                length() >0 )
#
#   expect_equal( create_gene_binary(mutation = mut_vt_na, snp_only = T) %>%
#                  ncol(), 0 )
#
#
# })
#
#
# # test include_silent arg----
# # add general tests
# # What happens  when Variant_Classification is NA for some samples in passed data? - Maybe need to add warning to tell user about NAs
# test_that("test include_silent arg", {
#
#   #general tests: input T or F (default is F)
#   expect_error( create_gene_binary(mutation=gnomeR::mut, include_silent = T), NA)
#
#   expect_warning( create_gene_binary(mutation=gnomeR::mut,
#                                      include_silent =  T, recode_aliases = FALSE), NA)
#
#
#   #What if NA for Variant_Classificaiton?
#   # note: without Variant Type, the create_gene_binary() still run without error
#   #       snp_only=F will provide full list results and snp_only=T will provide 0 col result
#
#   mut_vc_na<- gnomeR::mut %>%
#     dplyr::mutate(Variant_Classification=NA)
#
#   expect_equal( create_gene_binary(mutation = mut_vc_na, include_silent = F) %>% ncol(), 0 )
#
#   expect_true( create_gene_binary(mutation = mut_vc_na, include_silent = T) %>% ncol() > 0 )
#
# })
#
#
#
# # test rm_empty arg----
# # add general tests
# test_that("test rm_empty arg", {
#
#
#   #example test
#   expect_equal(TRUE, TRUE)
# })
#
#