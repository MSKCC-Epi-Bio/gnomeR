
# Test Binary Matrix Arguments -----------------------------------------------------------

# General tests ---
# test_that("test simplify marix with no errors", {
#
#   cna <- pivot_cna_longer(gnomeR::cna)
#   samples <- as.character(unique(gnomeR::mut$Tumor_Sample_Barcode))[1:50]
#   gen_dat <- create_gene_binary(samples = samples,
#                                 mutation = gnomeR::mut,
#                                 fusion = gnomeR::fusion,
#                                 cna = cna)
#
#   gen_dat2 <- gnomeR::summarize_by_patient(gen_dat)
#
#   tp53_1 <- gen_dat %>% select(contains("TP53"))
#   tp53_1 <- (tp53_1$TP53== 1 |tp53_1$TP53.Del== 1)%>% sum()
#
#   tp53_2 <- gen_dat2$TP53 %>% sum()
#   expect_equal(tp53_1, tp53_2)
#
#   expect_equal(length(table(gen_dat2$TP53)), 2)
# })

test_that("test simplify marix needs a data frame", {

  expect_error(summarize_by_patient(gene_binary = c(1:10)))

})

test_that("test simplify matrix needs a data frame", {

  samples <- as.character(unique(gnomeR::mutations$sampleId))[1:50]
  gen_dat <- create_gene_binary(samples = samples,
                                mutation = gnomeR::mutations,
                                fusion = gnomeR::sv)

  expect_no_error(gen_dat2 <- gnomeR::summarize_by_patient(gen_dat))


})

test_that("test output matches summarize_by_gene if only 1 sample per patient", {

  # restrict to 1 sample/patient
  samples <- gnomeR::mutations %>%
    mutate(patient_id = gnomeR::extract_patient_id(sampleId)) %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    pull(sampleId)

  gen_dat <- create_gene_binary(samples = samples,
                                mutation = gnomeR::mutations,
                                fusion = gnomeR::sv)

  expect_equal(gnomeR::summarize_by_gene(gen_dat) %>% select(-sample_id),
               gnomeR::summarize_by_patient(gen_dat) %>% select(-patient_id))


})

test_that("test that output is only 1 row/patient", {

  # restrict to 1 sample/patient
  # bind gnomeR::mutations to ensure there are duplicates
  samples <- bind_rows(gnomeR::mutations,
                       gnomeR::mutations) %>%
    mutate(patient_id = gnomeR::extract_patient_id(sampleId)) %>%
    pull(sampleId)

  gen_dat <- create_gene_binary(samples = samples,
                                mutation = gnomeR::mutations,
                                fusion = gnomeR::sv)

  n_rec_pt <- gnomeR::summarize_by_patient(gen_dat) %>%
    count(patient_id, name = "n_rec_pt") %>%
    distinct(n_rec_pt) %>%
    pull(n_rec_pt)

  expect_equal(n_rec_pt, n_rec_pt)
})

test_that("test that genes are properly summarized", {
  samples <- Reduce(intersect, list(gnomeR::mutations$sampleId, gnomeR::cna$sampleId,
                                    gnomeR::sv$sampleId))


  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                                   mutation = gnomeR::mutations,
                                                   cna = gnomeR::cna,
                                                   fusion = gnomeR::sv,
                                                   specify_panel = "impact") %>%
                   dplyr::select(c(sample_id, starts_with("ARI"), starts_with("MAPK1"), starts_with("ERG"))))


  sum_impact <- summarize_by_patient(bin_impact) %>%
    dplyr::mutate(across(!patient_id, as.numeric)) %>%
    arrange(patient_id)


  bin_impact_test <- bin_impact %>%
    mutate(patient_id = gnomeR::extract_patient_id(sample_id)) %>%
    select(-sample_id) %>%
    tidyr::pivot_longer(cols = c(everything(), -patient_id)) %>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x)) %>%
    mutate(across(.cols = c(-patient_id),
                  ~ifelse(. > 0, 1, 0))) %>%
    mutate_if(~ all(is.na(.)), ~as.numeric(NA_integer_))  %>%
    relocate(colnames(sum_impact)) %>%
    arrange(patient_id)

  expect_equal(sum_impact, bin_impact_test)
})

test_that("test alteration hierarchy", {
  # alteration if alteration on any sample
  # no alteration if no alteration on all samples (no NAs)
  # else NA

  # dummy binary matrix with comparisons of interest
  bin_impact_dummy <- tribble(~sample_id, ~ALK, ~EGFR, ~KRAS, ~TP53, ~KMT2D, ~SMARCA4,
                              paste0(gnomeR::mutations$patientId[1], "-T01", "-IM5"), 1, 0, 1, 1, 0, NA_integer_,
                              paste0(gnomeR::mutations$patientId[1], "-T01", "-IM5"), 1, 0, 0, NA_integer_, NA_integer_, NA_integer_)

  sum_impact <- summarize_by_patient(bin_impact_dummy) %>%
    dplyr::mutate(across(!patient_id, as.numeric)) %>%
    arrange(patient_id)

  bin_impact_expected <- tribble(~patient_id, ~ALK, ~EGFR, ~KRAS, ~TP53, ~KMT2D, ~SMARCA4,
                                 gnomeR::mutations$patientId[1], 1, 0, 1, 1, NA_integer_, NA_integer_) %>%
    mutate(across(-patient_id, as.numeric))

  expect_equal(sum_impact, bin_impact_expected)
})


test_that("test what happens to columns with all NA", {

  samples <- Reduce(intersect, list(gnomeR::mutations$sampleId, gnomeR::cna$sampleId,
                                    gnomeR::sv$sampleId))


  bin_impact <-  create_gene_binary(samples=samples,
                                    mutation = gnomeR::mutations,
                                    cna = gnomeR::cna,
                                    fusion = gnomeR::sv,
                                    specify_panel = "impact") %>%
    select(c(sample_id, starts_with("AR"), starts_with("PLCG2"), starts_with("PPM1D")))


  sum_impact <- summarize_by_patient(bin_impact)%>%
    mutate(across(!patient_id, as.numeric)) %>%
    arrange(patient_id)


  bin_impact_test <- bin_impact %>%
    mutate(patient_id = gnomeR::extract_patient_id(sample_id)) %>%
    select(-sample_id) %>%
    tidyr::pivot_longer(cols = c(everything(), -patient_id)) %>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x)) %>%
    mutate(across(.cols = c(-patient_id),
                  ~ifelse(. > 0, 1, 0))) %>%
    mutate_if(~ all(is.na(.)), ~as.numeric(NA_integer_))  %>%
    relocate(colnames(sum_impact)) %>%
    arrange(patient_id)

  expect_equal(sum_impact, bin_impact_test)

  expect_equal(ncol(as.data.frame(
    sum_impact[,sapply(sum_impact, function(x) all(is.na(x)))])), 1)

})


test_that("other vars are retained", {
  samples <- gnomeR::mutations %>%
    mutate(patient_id = gnomeR::extract_patient_id(sampleId)) %>%
    head(50) %>%
    pull(sampleId) %>%
    unique()

  bin_impact <- create_gene_binary(samples = samples,
                                   mutation = gnomeR::mutations,
                                   cna = gnomeR::cna,
                                   fusion = gnomeR::sv,
                                   specify_panel = "impact") %>%
    select(c(sample_id, starts_with("AR"), starts_with("PLCG2"), starts_with("PPM1D")))

  set.seed(20230828)

  bin_impact$random_color = sample(c("blue", "red", "yellow"),
                                   size = length(samples), replace = TRUE)

  expect_true("random_color" %in% names(bin_impact))
  sum_impact <- summarize_by_patient(bin_impact,
                                     other_vars = "random_color")

  expect_true("random_color" %in% names(sum_impact))
  expect_true("blue" %in% sum_impact$random_color)
})


test_that("no warning message thrown when only 1 alt type", {

  samples <- gnomeR::mutations$sampleId
  bin.mut <- create_gene_binary(
    samples = samples, mutation = gnomeR::mutations,
    mut_type = "omit_germline", snp_only = FALSE,
    include_silent = FALSE
  )

  expect_no_warning(summarize_by_patient(bin.mut))
})
