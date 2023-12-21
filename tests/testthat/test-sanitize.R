
# General Test ----------------------------------------------------------------
test_that("both sanitize functions run with no errors", {
  mut <- rename_columns(gnomeR::mutations)
  sv <- rename_columns(sv)

  expect_error(.sanitize_mutation_input(mut, include_silent = FALSE), NA)
  expect_error(.sanitize_mutation_input(mut, include_silent = TRUE), NA)
  expect_error(.sanitize_fusion_input(sv), NA)
})


test_that("test to see what happens if pass sanitize a vector", {
  mut <- rename_columns(gnomeR::mutations)
  sv <- rename_columns(sv)

  expect_error(.sanitize_mutation_input(gnomeR::mutations %>%
                                           select(-Hugo_Symbol),
                                       include_silent = FALSE))
  expect_error(.sanitize_fusion_input(gnomeR::sv %>%
                                         select(-Hugo_Symbol)))
})

test_that("alterations properly recoded using internal func", {
  cna <- gnomeR::cna[1:10,]%>%
    rename_columns() %>%
    .sanitize_cna_input()

  expect_true("amplification" %in% names(table(cna$alteration)))
  expect_true("deletion" %in% names(table(cna$alteration)))

  table <- as.data.frame(table(cna$alteration))

  expect_equal(table$Freq[table$Var1 == "deletion"], 7)
  expect_equal(table$Freq[table$Var1 == "amplification"], 3)
})


# --------------------------------------------------------------
# Test fusion in variant classification

test_that("test fusion in variant classification", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation$variant_classification[mutation$variant_classification == "In_Frame_Del"] <- "fusion"

  expect_error(.sanitize_mutation_input(mutation, include_silent = F), "It looks like you have fusions in your mutation data frame.*")
})


test_that("test fusion in variant classification", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-mutation_status)

  expect_warning(.sanitize_mutation_input(mutation, include_silent = F), "A mutation_status column*")

  mutation = mutation %>%
    mutate(mutation_status = "SOMATIC")

  expect_no_error(.sanitize_mutation_input(mutation, include_silent = F))
})

# variant type

test_that("test variant type", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-c(variant_type, reference_allele))

  expect_error(.sanitize_mutation_input(mutation, include_silent = F))
})

test_that("test variant type inference", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-c(variant_type))
  mutation$tumor_seq_allele_2 = mutation$reference_allele

  expect_warning(.sanitize_mutation_input(mutation, include_silent = F))
})


test_that("test variant type inference error", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation$tumor_seq_allele_2 = mutation$reference_allele
  mutation = mutation %>% select(-c(variant_type, reference_allele))

  expect_error(.sanitize_mutation_input(mutation, include_silent = F))
})


test_that("test that attributes are returned for input data col names", {

  mutation <- sanitize_mutation_input(gnomeR::mutations, include_silent = F)
  expect_true(!is.null(attr(mutation, "names_dict")))

  cna <- sanitize_cna_input(gnomeR::cna)
  expect_true(!is.null(attr(cna, "names_dict")))

  fus <- sanitize_fusion_input(gnomeR::sv)
  expect_true(!is.null(attr(fus, "names_dict")))


})



test_that("Check input colname is used in messaging", {

  mut_maf <- gnomeR::mutations %>%
    select(-variantType) %>%
    mutate(Tumor_Seq_Allele2 = NA)

  x <- capture_warning(sanitize_mutation_input(mut_maf, include_silent = FALSE))
  expect_true(str_detect(x$message, "Tumor_Seq_Allele2"))


})


# Test Specific Utils -----------------------------------------------------

test_that("clean and check columns works and cleans names", {

  expect_no_error(
    res <- .clean_and_check_cols(df_to_check = gnomeR::mutations))

  expect_true("hugo_symbol" %in% names(res))
})

test_that("clean and check columns returns character vectors", {

  mut <- gnomeR::mutations[1:10, ]
  mut$sampleId <- as.factor(mut$sampleId)

  expect_no_error(
    res <- .clean_and_check_cols(df_to_check = mut))

  expect_equal(class(res$sample_id), "character")

})

test_that("test check required cols works", {
  expect_no_error(.check_required_cols(rename_columns(gnomeR::mutations),
                       required_cols = c("sample_id", "hugo_symbol")))

})

test_that("test check required cols throws error when missing", {

  mut <- gnomeR::mutations %>%
    rename("sample_id" = sampleId) %>%
    select(-sample_id)

  expect_error(
    .check_required_cols(mut, required_cols = c("sample_id", "hugo_symbol"))
    )

})


