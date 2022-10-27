# # General -----------------------------------------------------------------
# test_that("works with basic input", {
#
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
#   gene_binary <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
#                           mut_type = "somatic_only", snp_only = FALSE) %>%
#     select(TP53, TP53.Del, APC, RB1, RB1.Del)
#
#   expect_error(tbl_genomic(gene_binary = gene_binary), NA)
# })
#
# test_that("pass both gene subset and freq", {
#
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
#   gene_binary <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
#                                  mut_type = "somatic_only", snp_only = FALSE) %>%
#     select(TP53, TP53.Del, APC, RB1, RB1.Del)
#
#   expect_message(tbl_genomic(gene_binary = gene_binary,
#                                    gene_subset = "TP53",
#                                    freq_cutoff = .1), "*")
# })
#
# test_that("test basic args", {
#
#   expect_error(tbl_genomic(gene_binary = c(1:10)))
#
#   samples <- mut$Tumor_Sample_Barcode[1:10]
#   gene_binary <- create_gene_binary(samples = samples,
#                                  mutation = mut,
#                                  cna = cna,
#                                  mut_type = "somatic_only",
#                                  snp_only = FALSE)  %>%
#     select(1:20)
#
#   expect_error(tbl_genomic(
#     gene_binary = gene_binary,
#     freq_cutoff = 2
#     ), "Please select a `freq_cutoff`")
#
#   # expect_error(tbl_genomic(
#   #   gene_binary = gene_binary,
#   #   by = c("TP53", "APC"),
#   #   freq_cutoff = 0
#   # ), "*")
#
#   expect_warning(tbl_genomic(
#     gene_binary = gene_binary,
#     gene_subset = c("TP53", "not_in_data")
#   ), "*")
#
#   expect_error(tbl_genomic(
#     gene_binary = gene_binary,
#     freq_cutoff = 1,
#     freq_cutoff_by_gene = FALSE
#   ), "No genes*")
#
#
# })
#
#
# # freq_cutoff -----------------------------------------------------------------
# test_that("check freq cutoff", {
#
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
#   gene_binary <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
#                                  mut_type = "somatic_only", snp_only = FALSE) %>%
#     select(TP53, TP53.Del, APC, RB1, RB1.Del)
#
#
#   # freq_cutoff_by_gene = FALSE
#   expect_error(by_gene <- tbl_genomic(gene_binary = gene_binary,
#                                    freq_cutoff = .2,
#                                    freq_cutoff_by_gene = FALSE), NA)
#
#   sums <- gene_binary %>%
#     tidyr::pivot_longer(everything()) %>%
#     group_by(name) %>%
#     summarise(sum_g = sum(value, na.rm = TRUE)/nrow(gene_binary))
#
#   over_cut<- sums %>% filter(sum_g >= .2) %>%
#     pull(name)
#
#   expect_equal(by_gene$table_body$variable, over_cut)
#
#   # freq_cutoff_by_gene = TRUE
#   expect_error(by_alt <- tbl_genomic(gene_binary = gene_binary,
#                                               freq_cutoff = .2,
#                                               freq_cutoff_by_gene = TRUE), NA)
#
#   over_cut <- sums %>%
#     tidyr::separate(col = name, into = c("name1", "name2"),
#                     remove =  FALSE, fill = "right") %>%
#     group_by(name1) %>%
#     mutate(s = sum(sum_g)) %>%
#     filter(s >= .2) %>%
#     pull(name) %>% unique()
#
#   expect_equal(sort(by_alt$table_body$variable), sort(over_cut))
#
# })
#
# # Test by variable in table ---------------------------------------------------
# test_that("test by variable not in data", {
#
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
#   gene_binary <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
#                                  mut_type = "somatic_only",
#                                  snp_only = FALSE) %>%
#     select(TP53, TP53.Del, APC, RB1, RB1.Del)
#
#   expect_error(tbl_genomic(gene_binary = gene_binary,
#                       freq_cutoff = .2,
#                       by = "nothing",
#                       freq_cutoff_by_gene = FALSE))
# })
#
# test_that("test by variable bare or string", {
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
#   gene_binary <- create_gene_binary(samples = samples, mutation = mut, cna = cna,
#                                  mut_type = "somatic_only",
#                                  snp_only = FALSE) %>%
#     select(TP53, TP53.Del, APC, RB1, RB1.Del)
#
#   gene_binary <- gene_binary %>%
#     mutate(sex = sample(x = c("M", "F"),
#                         size = nrow(gene_binary), replace = TRUE))
#
#   expect_error(t1 <- tbl_genomic(gene_binary = gene_binary,
#                       by = "sex",
#                       freq_cutoff = .2,
#                       freq_cutoff_by_gene = FALSE), NA)
#
#   expect_error(t2 <- tbl_genomic(gene_binary = gene_binary,
#                       by = sex,
#                       freq_cutoff = .2,
#                       freq_cutoff_by_gene = FALSE), NA)
#
#   expect_equal(t1, t2)
#   })
