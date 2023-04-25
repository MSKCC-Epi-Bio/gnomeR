
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

test_that("only accecpts tbl_gene_binary object", {
  fake <- data.frame(sample_id = c(rep("samp", 5)),
                     TERT = c(rep(1, 3), 0, NA))

  expect_error(summarize_by_gene(fake))

  binmat <- gnomeR::create_gene_binary(mutation = gnomeR::mutations[1:10,],
                                       cna = gnomeR::cna,
                                       fusion = gnomeR::sv[1:10,])

  expect_no_error(summarize_by_gene(binmat))

})

test_that("test simplify marix will not take other input", {

  expect_error(summarize_by_gene(c(1:10)))
  expect_error(summarize_by_gene(list(1:10)))
  expect_error(summarize_by_gene(c("test")))
  expect_error(summarize_by_gene(Sys.Date()))

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


  sum_impact <- summarize_by_gene(bin_impact)

  # manually
  bin_impact_test <- bin_impact %>%
    tidyr::pivot_longer(!sample_id)%>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x))%>%
    mutate(across(!sample_id, ~ifelse(. > 0, 1, 0)))%>%
    as.data.frame()%>%
    relocate(colnames(sum_impact))

  # set this manually to match expected classes
  class(bin_impact_test) <- c("tbl_gene_binary", class(bin_impact_test))

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


  sum_impact <- summarize_by_gene(bin_impact)


  bin_impact_test <- bin_impact %>%
    tidyr::pivot_longer(!sample_id)%>%
    mutate(name = str_remove(name, ".Amp|.Del|.fus"))%>%
    tidyr::pivot_wider(names_from = name, values_from = value, values_fn = function (x) sum(x))%>%
    mutate(across(!sample_id, ~ifelse(. > 0, 1, 0)))%>%
    as.data.frame()%>%
    relocate(colnames(sum_impact))%>%
    mutate_if(~ all(is.na(.)), ~as.numeric(NA_integer_))

  # set this manually to match expected classes
  class(bin_impact_test) <- c("tbl_gene_binary", class(bin_impact_test))

  expect_equal(sum_impact, bin_impact_test)

  x <- ncol(as.data.frame(
    sum_impact[,sapply(sum_impact, function(x) all(is.na(x)))]))

  expect_equal(x, 1)

  })



