
test_that("works with basic input", {

  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
  binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                          mut_type = "somatic_only", snp_only = FALSE) %>%
    select(TP53, TP53.Del, APC, RB1, RB1.Del)

  expect_error(genomic_tbl_summary(binary_matrix = binary_matrix), NA)
})


test_that("check freq cutoff", {

  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
  binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                                 mut_type = "somatic_only", snp_only = FALSE) %>%
    select(TP53, TP53.Del, APC, RB1, RB1.Del)

  sums <- binary_matrix %>%
    tidyr::pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(sum_g = sum(value, na.rm = TRUE)/nrow(binary_matrix))

  # HERE!!!!! THIS DOESN"T WORK
  expect_error(genomic_tbl_summary(binary_matrix = binary_matrix,
                                   freq_cutoff = .2, freq_cutoff_by_gene = FALSE), NA)

})
