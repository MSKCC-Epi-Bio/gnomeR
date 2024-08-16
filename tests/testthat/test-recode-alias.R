# test recode alias ----

test_that("missing columns of interest", {

  alias_table <- tibble::tribble(~hugo, ~alias,
                                 "CCND1",	"U21B31, BCL1, D11S287E, PRAD1")%>%
      dplyr::mutate(alias = as.list(strsplit(alias, ", ")))%>%
      tidyr::unnest(alias)

  genomic_df <- tibble::tribble(~hugo_symbol,
                                "CCND1",
                                "MYC",
                                "BCL1")

  expect_error(recode_alias(genomic_df, alias_table = alias_table), "Can't find required*")

})

test_that("warning if string of multiple aliases provided", {

  alias_table <- tibble::tribble(~hugo_symbol, ~alias,
                                 "CCND1",	"U21B31, BCL1, D11S287E, PRAD1",
                                 "RB1", "RB, PPP1R130, OSRC")
  alias_table2 <- tibble::tribble(~hugo_symbol, ~alias,
                                 "CCND1",	"U21B31,BCL1,D11S287E,PRAD1")

  alias_table3 <- tibble::tribble(~hugo_symbol, ~alias,
                                  "CCND1",	list("U21B31", "BCL1", "D11S287E", "PRAD1"))


  genomic_df <- tibble::tribble(~hugo_symbol,
                                "CCND1",
                                "MYC",
                                "BCL1")

  expect_error(recode_alias(genomic_df, alias_table), "Error with*")
  expect_error(recode_alias(genomic_df, alias_table2), "Error with*")
  expect_error(recode_alias(genomic_df, alias_table3), "Error with*")
})

test_that("aliases are recoded properly", {
  alias_table <- tibble::tribble(~hugo_symbol, ~alias,
                                 "CCND1",	"U21B31, BCL1, D11S287E, PRAD1")%>%
    dplyr::mutate(alias = as.list(strsplit(alias, ", "))) %>%
    tidyr::unnest(alias)

  genomic_df <- tibble::tribble(~hugo_symbol,
                                "CCND1",
                                "MYC",
                                "BCL1")

  expect_warning(genomic_df2 <- recode_alias(genomic_df, alias_table))
  expect_warning(genomic_df3 <- recode_alias(genomic_df))
  expect_warning(genomic_df4 <- recode_alias(genomic_df, "IMPACT"))

  test <- table(genomic_df2$hugo_symbol)
  test <- test[names(test) == "CCND1"]
  names(test) <- NULL

  expect_equal(class(test), "integer")
  expect_equal(test, 2)


  test <- table(genomic_df3$hugo_symbol)
  test <- test[names(test) == "CCND1"]
  names(test) <- NULL

  expect_equal(class(test), "integer")
  expect_equal(test, 2)

  test <- table(genomic_df4$hugo_symbol)
  test <- test[names(test) == "CCND1"]
  names(test) <- NULL

  expect_equal(class(test), "integer")
  expect_equal(test, 2)
})

test_that("aliases are recoded properly in create_gene_binary", {
  alias_table <- tibble::tribble(~hugo_symbol, ~alias,
                                 "CCND1",	"U21B31, BCL1, D11S287E, PRAD1")%>%
    dplyr::mutate(alias = as.list(strsplit(alias, ", "))) %>%
    tidyr::unnest(alias)


  mut <- head(gnomeR::mutations)
  mut$hugoGeneSymbol[1] <- "U21B31"

  cna <- head(gnomeR::cna)
  cna$hugoGeneSymbol[1] <- "BCL1"
  cna$hugoGeneSymbol[2] <- "U21B31"

  samples <- unique(c(mut$sampleId, cna$sampleId))
  expect_warning(genomic_df2 <- create_gene_binary(samples = samples,
                                                   mutation = mut,
                                                   cna = cna))
})
