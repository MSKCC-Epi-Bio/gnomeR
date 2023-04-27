# General -----------------------------------------------------------------
test_that("works with basic input", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                          mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")))

  expect_no_error(tbl_genomic(gene_binary = gene_binary))
})



# freq_cutoff -----------------------------------------------------------------
test_that("check freq cutoff", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only", snp_only = FALSE)


  # freq_cutoff_by_gene = FALSE
    by_gene <- gene_binary %>%
      subset_by_frequency(0.025) %>%
      tbl_genomic(gene_binary = .)

  sums <- gene_binary %>%
    select(-'sample_id') %>%
    tidyr::pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(sum_g = sum(value, na.rm = TRUE)/nrow(gene_binary))

  over_cut <- sums %>% filter(sum_g >= .025)

  vec <- over_cut$sum_g
  names(vec) <- over_cut$name

  vec <- sort(vec, decreasing = TRUE)

  expect_equal(sort(by_gene$table_body$variable), sort(names(vec)))


})

# Test by variable in table ---------------------------------------------------
test_that("test by variable not in data", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only",
                                 snp_only = FALSE) %>%
    select(all_of(c("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")))

  expect_error(tbl_genomic(gene_binary = gene_binary,
                      by = "nothing"))
})

test_that("test by variable bare or string", {
  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only",
                                 snp_only = FALSE) %>%
    select(all_of(c("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")))

  gene_binary <- gene_binary %>%
    mutate(sex = sample(x = c("M", "F"),
                        size = nrow(gene_binary), replace = TRUE))

  expect_error(t1 <- tbl_genomic(gene_binary = gene_binary,
                      by = "sex"), NA)

  expect_error(t2 <- tbl_genomic(gene_binary = gene_binary,
                      by = sex), NA)

  expect_equal(t1, t2)
  })


test_that("test ... to tbl_summary", {
  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only",
                                    snp_only = FALSE) %>%
    select(all_of(c("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")))

  expect_error(tbl_genomic(gene_binary = gene_binary,
                           statistic = list(everything() ~"{n}")), NA)


  })


test_that("you can pass gtsummary functions to tbl_genomic()",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))



  expect_no_error(tbl_genomic(gene_binary = gene_binary) %>%
                    gtsummary::bold_labels())

})


test_that("gene_binary is a tibble or dataframe",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))

  gene_binary = as.list(gene_binary)

  expect_error(tbl_genomic(gene_binary), "`gene_binary=` argument must be*")

})

test_that("check no duplicated IDs",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))

  gene_binary = rbind(gene_binary, gene_binary)

  expect_error(tbl_genomic(gene_binary), "Duplicate `sample_ids`*")

})

test_that("check freq_cutoff is not used",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))

  expect_error(tbl_genomic(gene_binary, freq_cutoff = 10), "The `freq_cutoff`*")

})

test_that("check freq_cutoff_by_gene is not used",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))

  expect_error(tbl_genomic(gene_binary, freq_cutoff_by_gene = TRUE), "The `freq_cutoff_by_gene`*")

})


test_that("check gene_subset is not used",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select(all_of(c('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')))

  expect_error(tbl_genomic(gene_binary, gene_subset = TRUE), "The `gene_subset`*")

})
