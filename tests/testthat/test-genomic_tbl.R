# General -----------------------------------------------------------------
test_that("works with basic input", {

  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
  binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                          mut_type = "somatic_only", snp_only = FALSE) %>%
    select(TP53, TP53.Del, APC, RB1, RB1.Del)

  expect_error(genomic_tbl_summary(binary_matrix = binary_matrix), NA)
})

test_that("pass both gene subset and freq", {

  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
  binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                                 mut_type = "somatic_only", snp_only = FALSE) %>%
    select(TP53, TP53.Del, APC, RB1, RB1.Del)

  expect_message(genomic_tbl_summary(binary_matrix = binary_matrix,
                                   gene_subset = "TP53",
                                   freq_cutoff = .1), "*")
})


# freq_cutoff -----------------------------------------------------------------
test_that("check freq cutoff", {

  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:10]
  binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
                                 mut_type = "somatic_only", snp_only = FALSE) %>%
    select(TP53, TP53.Del, APC, RB1, RB1.Del)


  # freq_cutoff_by_gene = FALSE
  expect_error(by_gene <- genomic_tbl_summary(binary_matrix = binary_matrix,
                                   freq_cutoff = .2,
                                   freq_cutoff_by_gene = FALSE), NA)

  sums <- binary_matrix %>%
    tidyr::pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(sum_g = sum(value, na.rm = TRUE)/nrow(binary_matrix))

  over_cut<- sums %>% filter(sum_g >= .2) %>%
    pull(name)

  expect_equal(by_gene$table_body$variable, over_cut)

  # freq_cutoff_by_gene = TRUE
  expect_error(by_alt <- genomic_tbl_summary(binary_matrix = binary_matrix,
                                              freq_cutoff = .2,
                                              freq_cutoff_by_gene = TRUE), NA)

  over_cut <- sums %>%
    tidyr::separate(col = name, into = c("name1", "name2"),
                    remove =  FALSE) %>%
    group_by(name1) %>%
    mutate(s = sum(sum_g)) %>%
    filter(s >= .2) %>%
    pull(name) %>% unique()

  expect_equal(sort(by_alt$table_body$variable), sort(over_cut))

})
