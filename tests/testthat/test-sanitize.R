
# General Test ----------------------------------------------------------------
test_that("both sanitize functions run with no errors", {
  expect_error(sanitize_mutation_input(gnomeR::mutations, include_silent = FALSE), NA)
  expect_error(sanitize_mutation_input(gnomeR::mutations, include_silent = TRUE), NA)
  expect_error(sanitize_fusion_input(gnomeR::sv), NA)
})


test_that("both sanitize functions run without warnings",{
  expect_warning(sanitize_mutation_input(gnomeR::mutations, include_silent = FALSE), NA)
  expect_warning(sanitize_mutation_input(gnomeR::mutations, include_silent = TRUE), NA)
  expect_warning(sanitize_fusion_input(gnomeR::sv), NA)
})

test_that("test to see what happens if pass sanitize a vector", {
  expect_error(sanitize_mutation_input(gnomeR::mutations %>%
                                           select(-Hugo_Symbol),
                                       include_silent = FALSE))
  expect_error(sanitize_fusion_input(gnomeR::sv %>%
                                         select(-Hugo_Symbol)))
})

test_that("alterations properly recoded using internal func", {
  cna <- gnomeR::cna[1:10,]%>%
    sanitize_cna_input()

  expect_true("amplification" %in% names(table(cna$alteration)))
  expect_true("deletion" %in% names(table(cna$alteration)))

  table <- as.data.frame(table(cna$alteration))

  expect_equal(table$Freq[table$Var1 == "deletion"], 7)
  expect_equal(table$Freq[table$Var1 == "amplification"], 3)
})

# Test Required Columns -------------------------------------------------------
test_that("both sanitize functions run with no errors", {

  mutations <- select(gnomeR::mutations, -"hugoGeneSymbol")
  expect_error(sanitize_mutation_input(mutations), "The following*")

})

# --------------------------------------------------------------
## added by cw on 4/26/23
# sanitize_fusion_input(fusion)
# colnames(fusion)

# test fusion in variant classification

test_that("test fusion in variant classification", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation$variant_classification[mutation$variant_classification == "In_Frame_Del"] <- "fusion"

  expect_error(sanitize_mutation_input(mutation, include_silent = F), "It looks like you have fusions in your mutation data frame.*")
})

# check suggested columns
# mutation status column

test_that("test fusion in variant classification", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-mutation_status)

  expect_warning(sanitize_mutation_input(mutation, include_silent = F), "A mutation_status column*")

  mutation = mutation %>%
    mutate(mutation_status = "SOMATIC")

  expect_no_error(sanitize_mutation_input(mutation, include_silent = F))
})

# variant type

test_that("test variant type", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-c(variant_type, reference_allele))

  expect_error(sanitize_mutation_input(mutation, include_silent = F))
})

test_that("test variant type inference", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation = mutation %>% select(-c(variant_type))
  mutation$tumor_seq_allele2 = mutation$reference_allele

  expect_warning(sanitize_mutation_input(mutation, include_silent = F))
})


test_that("test variant type inference error", {
  mutation = gnomeR::mutations
  mutation <- rename_columns(mutation)
  column_names <- colnames(mutation)
  mutation$tumor_seq_allele2 = mutation$reference_allele
  mutation = mutation %>% select(-c(variant_type, reference_allele))

  expect_error(sanitize_mutation_input(mutation, include_silent = F))
})



