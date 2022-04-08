
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
# NOTE - SEE EDIT NEEDED IN ISSUE 150 before testing: https://github.com/MSKCC-Epi-Bio/gnomeR/issues/150
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
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
test_that("test include_silent arg", {

  #general tests: input T or F (default is F)
  expect_error( binary_matrix(mutation=gnomeR::mut, include_silent = T), NA)

  expect_warning( binary_matrix(mutation=gnomeR::mut, include_silent =  T), NA)


  #What if NA for Variant_Classificaiton?
  # note: without Variant Type, the binary_matrix() still run without error
  #       snp_only=F will provide full list results and snp_only=T will provide 0 col result

  mut_vc_na<- gnomeR::mut %>%
    dplyr::mutate(Variant_Classification=NA)

  expect_equal( binary_matrix(mutation = mut_vc_na, include_silent = F) %>% ncol(), 0 )

  expect_true( binary_matrix(mutation = mut_vt_na, include_silent = T) %>% ncol() > 0 )

})

# test cna_binary arg----
# add general tests
# I don't have an example of data that has cna values that aren't just 1 or 2. It would be helpful to
# find an example of data to test this using the API {cbioportalR}. Then make it smaller (just a few rows) and test using that
test_that("test for cna_binary arg", {

  # add general tests (default is T)
      ## If T, then the output should be all 0 or 1
  expect_identical(binary_matrix(cna = gnomeR::cna) %>%
                     purrr::map_dbl(~any(!(.x %in% c(0,1))) ) %>%
                      sum(), 0)
     ## If F, each column will represent only one different gene
  res_cna<- names(binary_matrix(cna = cna, cna_binary = F)) %>%
             stringr::str_replace(c(".cna"),"")

  expect_equal(length( unique(res_cna)),  length(res_cna))

  #example test
  expect_equal(TRUE, TRUE)
})

# test cna_relax arg----
# add general tests
# find an example of data to test this using the API {cbioportalR} that has both 1 and 2 values. Then make it smaller (just a few rows) and test using that
test_that("test for cna_relax arg", {

  # add general tests (default is F)
  expect_error(binary_matrix(cna=cna, cna_relax = T), NA)

  # Use a fake data to test if T then consider both 1 and -1 as 2 and -2
  cna_fake <- data.frame(gnomeR::cna[1:5,1],matrix(sample(seq(-2,2),5*20, replace=TRUE),nrow=5))
  names(cna_fake)<-names(gnomeR::cna)[1:21]

  amp.del.ct<-function(input_vec, amp_val, del_val){
    amp_ct<-sapply(input_vec, function(x){ as.numeric( x %in% amp_val )}) %>% sum()
    del_ct<-sapply(input_vec, function(x){ as.numeric( x %in% del_val )}) %>% sum()
    return(c(amp_ct, del_ct))
  }

  expect_equal( binary_matrix(cna=cna_fake, cna_relax=T) %>%
                  sapply(sum) %>%
                    as.vector(),
                apply(cna_fake[,-1], 1, amp.del.ct, amp_val=c(1,2), del_val=c(-1,-2) ) %>%
                        as.vector() )

   ### Note: the cna_relax=F case not working for ".Del"
  # expect_equal( binary_matrix(cna=cna_fake, cna_relax=F) %>%
  #                 sapply(sum) %>%
  #                 as.vector(),
  #               apply(cna_fake[,-1], 1, amp.del.ct, amp_val=c(2), del_val=c(-2) ) %>%
  #                 as.vector() )



  #example test
  expect_equal(TRUE, TRUE)
})



# test rm_empty arg----
# add general tests
test_that("test rm_empty arg", {


  #example test
  expect_equal(TRUE, TRUE)
})


