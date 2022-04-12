
# Test Binary Matrix Arguments -----------------------------------------------------------

# General tests ---
test_that("test binary_matrix with mutation runs with no errors", {

  expect_error(binary_matrix(mutation = gnomeR::mut), NA)

  res_mut <- binary_matrix(mutation = gnomeR::mut)
  expect_true(length(res_mut) > 0)
})


test_that("test binary_matrix with cna runs with no errors", {

  expect_error(binary_matrix(cna = gnomeR::cna), NA)

  res <- binary_matrix(cna = gnomeR::cna)
  expect_true(length(res) > 0)
})


test_that("test binary_matrix with fusions runs with no errors", {

  expect_error(binary_matrix(fusion = gnomeR::fusion), NA)

  res_fusion <- binary_matrix(fusion = gnomeR::fusion)
  expect_true(length(res_fusion) > 0)

  length(unique(gnomeR::fusion))
})

test_that("check cna with no alterations are omitted from results", {

  res <- binary_matrix(mutation = gnomeR::mut,
                       cna = gnomeR::cna,
                       fusion = gnomeR::fusion)
  cna_ids <- names(gnomeR::cna)[-1] %>%
    str_replace_all(fixed("."), "-")

  omitted_ids <- setdiff(cna_ids, rownames(res))

  omitted_ids <- omitted_ids %>%
    str_replace_all(fixed("-"), fixed("."))

  check_they_are_zero <- gnomeR::cna %>% select(all_of(omitted_ids)) %>%
    purrr::map_dbl(., ~sum(.x))

  expect_true(sum(check_they_are_zero) == 0)
})

# test samples argument ----
# what happens when you pass a vector? What about if you don't specify it (don't pass anything)?
# what happens when you pass impact samples (-IM5/IM6/IM7 etc)?  non impact samples? A mix?


test_that("Check binary_matrix() provide specific sample data if pass a vector", {

  #what happens when you pass a vector?
  mut_valid_sample_ids<-binary_matrix( mutation= gnomeR::mut) %>%
                     rownames() %>%
                     head(n=10)

  expect_equal(
    binary_matrix(sample=mut_valid_sample_ids, mutation=gnomeR::mut) %>%
      nrow(),
    length(mut_valid_sample_ids))

  #what about if you don't specify it (don't pass anything)?

  expect_lte(
    binary_matrix(mutation=gnomeR::mut) %>%
      nrow(),
    length(gnomeR::mut[['Tumor_Sample_Barcode']]))

})

test_that("Check binary_matrix() if sample entered in `sampl_id` with zero mutations/fusions/cna", {

  #what happens when you pass a vector?
  mut_valid_sample_ids<-binary_matrix( mutation= gnomeR::mut) %>%
    rownames() %>%
    head(n=10)

  add_no_mut_sample <- c(mut_valid_sample_ids, "no_mutations_fake_sample")
  binary_matrix_with_zero <-  binary_matrix(sample=add_no_mut_sample, mutation=gnomeR::mut)

  sum(binary_matrix_with_zero[nrow(binary_matrix_with_zero), ])
  expect_equal(
    sum(binary_matrix_with_zero[nrow(binary_matrix_with_zero), ]), 0)

  # should be one more obs in data frame with samples arg specified
  expect_equal(nrow(binary_matrix_with_zero) -1, length(mut_valid_sample_ids))

  # with no fusions in select sample---------
  binary_matrix_with_zero <-  binary_matrix(samples=add_no_mut_sample,
                                            mutation = gnomeR::mut,
                                            cna = gnomeR::cna,
                                            fusion = gnomeR::fusion)

  expect_false(any(str_detect(names(binary_matrix_with_zero), ".fus")))
  expect_equal(nrow(binary_matrix_with_zero), length(add_no_mut_sample))


  # with no cna in select sample---------
  cna_samp <- cna[, c(1, 100)]
  binary_matrix_with_zero <-  binary_matrix(samples=add_no_mut_sample,
                                            mutation = gnomeR::mut,
                                            cna = cna_samp,
                                            fusion = gnomeR::fusion)
  expect_false(any(str_detect(names(binary_matrix_with_zero), ".Amp")))
  expect_false(any(str_detect(names(binary_matrix_with_zero), "Del")))
  expect_false(any(str_detect(names(binary_matrix_with_zero), ".cna")))
  expect_equal(nrow(binary_matrix_with_zero), length(add_no_mut_sample))

})

# test with and without mut/fusion/cna args passed ----
# Functions should work with any one of the three passed
# Does it return results as expected?
# Trying with and without passing samples arg as well- does it return what you'd expect?
# what if mut is passed but doesn't have any rows? no columns? no rows or cols?
test_that("test inputting mut/fusion/cna args can leads to a data.frame output", {

  #Can we obtaine correct result format when either mut/fusion/cna inputted
  expect_error(binary_matrix())

  expect_true( binary_matrix( mutation = gnomeR::mut) %>% is.data.frame())

  expect_true( binary_matrix( fusion = gnomeR::fusion) %>% is.data.frame())

  expect_true( binary_matrix( cna = gnomeR::cna) %>% is.data.frame())

  # What if there is no row/col in the file passing to mutation

  #note: there is no error message if a inputting mut data is 0 rows;
  #      it only return a 0 row and 0 column result
  expect_equal( gnomeR::mut %>%
                  dplyr::filter(.data$Hugo_Symbol=="XXXXXXXXX") %>%
                  binary_matrix(mutation = .) %>%
                  nrow(), 0)

  expect_error( gnomeR::mut %>%
                dplyr::select(NULL) %>%
                binary_matrix(mutation = .))


})

# test mut_type argument ----

test_that("test incorrectly specified arg", {

  expect_error(binary_matrix(mutation = mut2,mut_type = "somatic_only",
                             specify_panel = "no"))
})


test_that("test inclusion of NAs in mut_type ", {
  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""

  #example test
  expect_warning(binary_matrix(mutation = mut2, specify_panel = "no"), "15 mutations*")
})



test_that("test inclusion of NAs in mut_type ", {

  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""

  # NA included by default (germline_omitted)
  expect_warning(see <- binary_matrix(mutation = mut2, specify_panel = "no"))
  check <-see$TP53[which(rownames(see)=="P-0000062-T01-IM3")]
  expect_equal(check, 1)

})

test_that("test inclusion of NAs in mut_type ", {
  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""


  # NA included with all
  see = binary_matrix(mutation = mut2, specify_panel = "no", mut_type = "all")
  expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],1)


  # NA no longer included with somatic_only
  see = binary_matrix(mutation = mut2, mut_type = "somatic_only", specify_panel = "no")
  expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],0)

  # NA no longer included with germline_only
  see = binary_matrix(mutation = mut2, mut_type = "germline_only", specify_panel = "no")
  expect_equal(ncol(see), 0)


})

# test snp_only arg----
# add general tests
# What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs
test_that("test the snp_only arg", {

  #general tests: input T or F (default is F)
  expect_error( binary_matrix(mutation=gnomeR::mut, snp_only = T), NA)

  expect_warning( binary_matrix(mutation=gnomeR::mut, snp_only = T), NA)

  #0 col return if 0 SNP been inputted or there is no SNP category in Variat Type been inputted
  mut_snp_zero<- gnomeR::mut %>%
                 dplyr::filter(.data$Variant_Type!="SNP")

  mut_snp_na<-mut_snp_zero
  mut_snp_na$Variant_Type<- droplevels(mut_snp_na$Variant_Type)

  expect_equal( binary_matrix(mutation=mut_snp_zero, snp_only = T) %>%
                ncol(),
                0)

  expect_equal( binary_matrix(mutation=mut_snp_na, snp_only = T) %>%
                  ncol(),
                0)

  #What if NA for Variant Type?
  # note: without Variant Type, the binary_matrix() still run without error
  #       snp_only=F will provide full list results and snp_only=T will provide 0 col result

  mut_vt_na<- gnomeR::mut %>%
              dplyr::mutate(Variant_Type=NA)

  expect_true( binary_matrix(mutation = mut_vt_na, snp_only = F) %>%
               length() >0 )

  expect_equal( binary_matrix(mutation = mut_vt_na, snp_only = T) %>%
                 ncol(), 0 )


})


# test include_silent arg----
# add general tests
# What happens  when Variant_Classification is NA for some samples in passed data? - Maybe need to add warning to tell user about NAs
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test cna_binary arg----
# add general tests
# I don't have an example of data that has cna values that aren't just 1 or 2. It would be helpful to
# find an example of data to test this using the API {cbioportalR}. Then make it smaller (just a few rows) and test using that
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test cna_relax arg----
# add general tests
# find an example of data to test this using the API {cbioportalR} that has both 1 and 2 values. Then make it smaller (just a few rows) and test using that
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})



# test rm_empty arg----
# add general tests
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})


