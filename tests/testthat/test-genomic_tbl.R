# General -----------------------------------------------------------------
test_that("works with basic input", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                          mut_type = "somatic_only", snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  expect_no_error(tbl_genomic(gene_binary = gene_binary, freq_cutoff = 0))
})

test_that("pass both gene subset and freq", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples,
                                    mutation = gnomeR::mutations, cna = gnomeR::cna,
                                 mut_type = "somatic_only", snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  expect_message(tbl_genomic(gene_binary = gene_binary,
                                   gene_subset = "MYC",
                                   freq_cutoff = .1), "*")


  expect_no_error(tbl_genomic(gene_binary = gene_binary,
                             gene_subset = "MYC",
                             freq_cutoff = .1))
})

test_that("test basic args", {

  expect_error(tbl_genomic(gene_binary = c(1:10)))

  samples <- mutations$sampleId[1:10]
  gene_binary <- create_gene_binary(samples = samples,
                                 mutation = mutations,
                                 cna = cna,
                                 mut_type = "somatic_only",
                                 snp_only = FALSE)  %>%
    select(1:20)

  expect_error(tbl_genomic(
    gene_binary = gene_binary,
    freq_cutoff = 2
    ), "Please select a `freq_cutoff`")

  expect_error(tbl_genomic(
    gene_binary = gene_binary,
    by = c("TP53", "APC"),
    freq_cutoff = 0
  ), "*")

  #broken here
  # expect_error(tbl_genomic(
  #   gene_binary = gene_binary,
  #   gene_subset = c("TP53", "not_in_data")
  # ), "*")

  expect_error(tbl_genomic(
    gene_binary = gene_binary,
    freq_cutoff = 1,
    freq_cutoff_by_gene = FALSE
  ), "No genes*")


})


# freq_cutoff -----------------------------------------------------------------
test_that("check freq cutoff", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only", snp_only = FALSE)  %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")


  # freq_cutoff_by_gene = FALSE
  expect_error(by_gene <- tbl_genomic(gene_binary = gene_binary,
                                   freq_cutoff = .025,
                                   freq_cutoff_by_gene = FALSE), NA)

  sums <- gene_binary %>%
    select(-'sample_id') %>%
    tidyr::pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(sum_g = sum(value, na.rm = TRUE)/nrow(gene_binary))

  over_cut<- sums %>% filter(sum_g >= .025) %>%
    pull(name)

  expect_equal(by_gene$table_body$variable, over_cut)

  # freq_cutoff_by_gene = TRUE
  expect_error(by_alt <- tbl_genomic(gene_binary = gene_binary,
                                              freq_cutoff = .025,
                                              freq_cutoff_by_gene = TRUE), NA)

  over_cut <- sums %>%
    tidyr::separate(col = name, into = c("name1", "name2"),
                    remove =  FALSE, fill = "right") %>%
    group_by(name1) %>%
    mutate(s = sum(sum_g)) %>%
    filter(s >= .025) %>%
    pull(name) %>% unique()

  expect_equal(sort(by_alt$table_body$variable), sort(over_cut))

})

# Test by variable in table ---------------------------------------------------
test_that("test by variable not in data", {

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only",
                                 snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  expect_error(tbl_genomic(gene_binary = gene_binary,
                      freq_cutoff = .025,
                      by = "nothing",
                      freq_cutoff_by_gene = FALSE))
})

test_that("test by variable bare or string", {
  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                 mut_type = "somatic_only",
                                 snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  gene_binary <- gene_binary %>%
    mutate(sex = sample(x = c("M", "F"),
                        size = nrow(gene_binary), replace = TRUE))

  expect_error(t1 <- tbl_genomic(gene_binary = gene_binary,
                      by = "sex",
                      freq_cutoff = .025,
                      freq_cutoff_by_gene = FALSE), NA)

  expect_error(t2 <- tbl_genomic(gene_binary = gene_binary,
                      by = sex,
                      freq_cutoff = .025,
                      freq_cutoff_by_gene = FALSE), NA)

  expect_equal(t1, t2)
  })



test_that("test ... to tbl_summary", {
  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only",
                                    snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  expect_error(tbl_genomic(gene_binary = gene_binary,
                           statistic = list(all_categorical() ~"{n}")), "*")


  })

test_that("you need to load gtsummary for ...",{

  suppressWarnings(library(gtsummary))
  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only",
                                    snp_only = FALSE) %>%
    select("FGFR1.Amp", "SOX17.Amp", "MYC.Amp", "MYC", "sample_id")

  #need to load package for `...`
  expect_no_error(tbl_genomic(gene_binary = gene_binary,
                              freq_cutoff = .025,
                              statistic = list(all_categorical() ~"{n}")))

})

test_that("you can pass gtsummary functions to tbl_genomic()",{

  samples <- as.character(unique(mutations$sampleId))[1:10]
  gene_binary <- create_gene_binary(samples = samples, mutation = mutations, cna = cna,
                                    mut_type = "somatic_only", snp_only = FALSE) %>%
    select('SMAD2', 'FGFR1.Amp', 'AKT1', 'SOX17.Amp', 'MYC', 'MYC.Amp', 'sample_id')



  expect_no_error(tbl_genomic(gene_binary = gene_binary, freq_cutoff = 0.025) %>%
                    gtsummary::bold_labels())

})

