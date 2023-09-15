
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
#   gen_dat2 <- gnomeR::summarize_by_gene(gen_dat)
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

  expect_error(summarize_by_gene(gene_binary = c(1:10)))

})

test_that("test simplify matrix needs a data frame", {

    samples <- as.character(unique(gnomeR::mutations$sampleId))[1:50]
    gen_dat <- create_gene_binary(samples = samples,
                                  mutation = gnomeR::mutations,
                                  fusion = gnomeR::sv)

    expect_no_error(gen_dat2 <- gnomeR::summarize_by_gene(gen_dat))


})


test_that("test that genes are properly summarized", {
  samples <- Reduce(intersect, list(gnomeR::mutations$sampleId, gnomeR::cna$sampleId,
                                    gnomeR::sv$sampleId))


  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                    mutation = gnomeR::mutations,
                                    cna = gnomeR::cna,
                                    fusion = gnomeR::sv,
                                    specify_panel = "impact")%>%
    dplyr::select(c(sample_id, starts_with("ARI"), starts_with("MAPK1"), starts_with("ERG"))))


  sum_impact <- summarize_by_gene(bin_impact)%>%
    dplyr::mutate(across(!sample_id, as.numeric))


  bin_impact_test <- bin_impact %>%
    tidyr::pivot_longer(!sample_id)%>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x))%>%
    mutate(across(!sample_id, ~ifelse(. > 0, 1, 0)))%>%
    relocate(colnames(sum_impact))

  expect_equal(sum_impact, bin_impact_test)


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


  sum_impact <- summarize_by_gene(bin_impact)%>%
    mutate(across(!sample_id, as.numeric))


  bin_impact_test <- bin_impact %>%
    tidyr::pivot_longer(!sample_id)%>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x))%>%
    mutate(across(!sample_id, ~ifelse(. > 0, 1, 0)))%>%
    relocate(colnames(sum_impact))%>%
    mutate_if(~ all(is.na(.)), ~as.numeric(NA_integer_))

  expect_equal(sum_impact, bin_impact_test)

  expect_equal(ncol(as.data.frame(
    sum_impact[,sapply(sum_impact, function(x) all(is.na(x)))])), 1)

  })


# test_that("all columns must be numeric to continue", {
#
# })

test_that("other vars are retained", {
  samples <- Reduce(intersect, list(gnomeR::mutations$sampleId,
                                    gnomeR::cna$sampleId,
                                    gnomeR::sv$sampleId))


  bin_impact <- create_gene_binary(samples = samples,
                                    mutation = gnomeR::mutations,
                                    cna = gnomeR::cna,
                                    fusion = gnomeR::sv,
                                    specify_panel = "impact") %>%
    select(c(sample_id, starts_with("AR"), starts_with("PLCG2"), starts_with("PPM1D")))

  set.seed(20230828)

  bin_impact$random_color = sample(c("blue", "red", "yellow"),
                                   size = 50, replace = TRUE)

  expect_true("random_color" %in% names(bin_impact))
  sum_impact <- summarize_by_gene(bin_impact,
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

  expect_no_warning(summarize_by_gene(bin.mut))



})





