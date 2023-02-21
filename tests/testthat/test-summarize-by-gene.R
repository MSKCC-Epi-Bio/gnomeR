
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

test_that("test simplify marix needs a data frame", {

    samples <- as.character(unique(gnomeR::mutations$sampleId))[1:50]
    gen_dat <- create_gene_binary(samples = samples,
                                  mutation = gnomeR::mutations,
                                  fusion = gnomeR::sv)

    expect_no_error(gen_dat2 <- gnomeR::summarize_by_gene(gen_dat))


})




