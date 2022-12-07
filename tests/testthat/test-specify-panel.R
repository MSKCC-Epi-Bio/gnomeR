# test specify_impact_panels() ----

test_that("gene binary with muts with impact specified", {
  samples <- unique(gnomeR::mutations$sampleId)[1:10]

  test_gene <- gnomeR::mutations %>%
    filter(sampleId %in% samples) %>%
    mutate(hugoGeneSymbol = "XXXTEST")%>%
    group_by(sampleId)%>%
    filter(row_number()==1) #select only 1 of each sample

  mut_test <- gnomeR::mutations %>%
    filter(sampleId %in% samples)%>%
    rbind(test_gene)

  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                           mutation = mut_test,
                                           specify_panel = "impact"),
                 "1 column*")


  num_genes <- length(unique(c(mut_test$hugoGeneSymbol)))

  expect_true(ncol(bin_impact) == 1 + num_genes)

  expect_true(nrow(bin_impact) == length(samples))

  expect_true("XXXTEST" %in% colnames(bin_impact))

  expect_equal(class(bin_impact$XXXTEST), "logical")

  testNA <- table(is.na(bin_impact$XXXTEST)) %>%
    as.data.frame()%>%
    filter(Var1 == "TRUE")

  expect_equal(nrow(bin_impact), testNA$Freq)
})

test_that("gene binary test cna with impact specified", {
  samples <- unique(gnomeR::cna$sampleId)[1:10]

  bin_impact <-  create_gene_binary(samples=samples,
                                   cna = gnomeR::cna,
                                   specify_panel = "impact")


  cna_test <- gnomeR::cna %>%
    filter(sampleId %in% samples)

  num_genes <- length(unique(c(cna_test$hugoGeneSymbol)))



  expect_true(ncol(bin_impact) == 1 + num_genes)

  expect_true(nrow(bin_impact) == length(samples))


})

test_that("gene binary test fusion with impact specified", {
  samples <- unique(gnomeR::sv$sampleId)[1:10]

  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                   fusion = gnomeR::sv,
                                   specify_panel = "impact"),
                 "* have all missing*")


  sv_test <- gnomeR::sv %>%
    filter(sampleId %in% samples)

  #some fusions have NA listed in second site so drop from list
  num_genes <- unique((c(sv_test$site1HugoSymbol, sv_test$site2HugoSymbol)))%>%
    na.omit()%>%
    length()

  expect_true(ncol(bin_impact) == 1 + num_genes)

  expect_true(nrow(bin_impact) == length(samples))

  })

test_that("gene binary with all three types of alt and impact only",{
  #select samples that are in all three types of mutations
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                       gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- unique(samples)[1:5]

  bin_impact <-  create_gene_binary(samples=samples,
                                   mutation = gnomeR::mutations,
                                   cna = gnomeR::cna,
                                   fusion = gnomeR::sv,
                                   specify_panel = "impact")

  mut_test <- gnomeR::mutations %>%
    filter(sampleId %in% samples)

  cna_test <- gnomeR::cna %>%
    filter(sampleId %in% samples)

  sv_test <- gnomeR::sv %>%
    filter(sampleId %in% samples)

  mut_genes <- length(unique(mut_test$hugoGeneSymbol))
  cna_genes <- length(unique(cna_test$hugoGeneSymbol))
  sv_genes <- unique(c(sv_test$site1HugoSymbol,
                              sv_test$site2HugoSymbol))%>%
    na.omit() %>%
    length()

  expect_true(ncol(bin_impact) == 1 + mut_genes + cna_genes + sv_genes)

  expect_true(nrow(bin_impact) == length(samples))


})

test_that("test 0 impact genes", {


  mut_test <- gnomeR::mutations %>%
    mutate(sampleId = patientId) #removes IM IH endings

  cna_test <- gnomeR::cna %>%
    mutate(sampleId = patientId)

  sv_test <- gnomeR::sv %>%
    mutate(sampleId = patientId)


  samples <- mut_test$sampleId[mut_test$sampleId %in%
                                          cna_test$sampleId]
  samples <-  samples[samples %in% sv_test$sampleId]
  samples <- unique(samples)[1:5]

  expect_error(bin_impact <-  create_gene_binary(samples=samples,
                                   mutation = mut_test,
                                   cna = cna_test,
                                   fusion = sv_test,
                                   specify_panel = "impact"),
               "There are no IMPACT*")



})

test_that("endings don't match impact list", {
  #removes true IM IH endings and adds IM11 which doesn't exist

  cna_test <- gnomeR::cna %>%
    mutate(sampleId = paste0(patientId, "-T01-IM2"))%>%
    head()

  sv_test <- gnomeR::sv %>%
    mutate(sampleId = paste0(patientId, "-T01-IMA"))%>%
    head()


  patients <- unique(gnomeR::mutations$patientId)[1:10] #pick patient here because will mutate later
  mut_test1 <- gnomeR::mutations %>%
    filter(patientId %in% patients) %>%
    mutate(hugoGeneSymbol = "XXXTEST")%>%
    group_by(sampleId)%>%
    filter(row_number()==1)%>%
    mutate(sampleId = paste0(patientId, "-T01-IM11"))

  mut_test <- gnomeR::mutations%>%
    filter(patientId %in% patients)%>%
    mutate(sampleId = paste0(patientId, "-T01-IM11"))%>%
    rbind(mut_test1)

  samples <- unique(c(mut_test$sampleId, cna_test$sampleId,
                    sv_test$sampleId))

  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                                 mutation = mut_test,
                                                 cna = cna_test,
                                                 fusion = sv_test,
                                                 specify_panel = "impact"),
                 "Couldn't infer IMPACT*")


  #table should still be created just no NAs filled in for non-IMPACT genes
  mut_genes <- length(unique(mut_test$hugoGeneSymbol))
  cna_genes <- length(unique(cna_test$hugoGeneSymbol))
  sv_genes <- unique(c(sv_test$site1HugoSymbol,
                       sv_test$site2HugoSymbol))%>%
    na.omit() %>%
    length()

  expect_true(ncol(bin_impact) == 1 + mut_genes + cna_genes + sv_genes)

  expect_true(nrow(bin_impact) == length(samples))

  summary <- sum(bin_impact$XXXTEST)

  expect_equal(summary, length(patients))
})


test_that("check df for specify panel has correct names",{
  #select samples that are in all three types of mutations
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                                          gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- unique(samples)[1:5]

  sp <- as.data.frame(samples) %>%
    mutate(panels = "IMPACT341")

  mut <- gnomeR::mutations %>%
    filter(sampleId %in% samples)

  expect_error(bin_impact <-  create_gene_binary(samples=samples,
                                    mutation = mut,
                                    cna = gnomeR::cna,
                                    fusion = gnomeR::sv,
                                    specify_panel = sp), "Dataframe*")


})


test_that("check specify_panel with dataframe works",{
  #select samples that are in all three types of mutations
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                                          gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- unique(samples)[1:5]

  sp <- as.data.frame(samples) %>%
    mutate(panels = "IMPACT341")

  mut <- gnomeR::mutations %>%
    filter(sampleId %in% samples)
  mut[1, 1] <- "ZZZZ"

  sp <- rename(sp,
               "sample_id" = samples,
               "panel_id" = panels)

  expect_message(bin_impact <-  create_gene_binary(samples=samples,
                                                   mutation = mut,
                                                   cna = gnomeR::cna,
                                                   fusion = gnomeR::sv,
                                                   specify_panel = sp), "*")

  # recode so one is on later panel with CSF3R

  sp[1, 2] <- "IMPACT468"
  expect_message(bin_impact2 <-  create_gene_binary(samples=samples,
                                                 mutation = mut,
                                                 cna = gnomeR::cna,
                                                 fusion = gnomeR::sv,
                                                 specify_panel = sp), "*")

  expect_gt(sum(is.na(bin_impact$CSF3R)), sum(is.na(bin_impact2$CSF3R)))


})


test_that("check specify_panel with dataframe that doesn't have all samples",{
  #select samples that are in all three types of mutations
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                                          gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- unique(samples)[1:5]

  sp <- as.data.frame(samples) %>%
    mutate(panels = "IMPACT341")

  mut <- gnomeR::mutations %>%
    filter(sampleId %in% samples)


  sp <- rename(sp,
               "sample_id" = samples,
               "panel_id" = panels)

  sp <- head(sp, 1)

  expect_no_message(bin_impact <-  create_gene_binary(samples=samples,
                                                   mutation = mut,
                                                   cna = gnomeR::cna,
                                                   fusion = gnomeR::sv,
                                                   specify_panel = sp))

  expect_equal(sum(apply(bin_impact, 1, function(x) sum(is.na(x))) > 0), 1)



})


test_that("Returns error if missing in specify panel DF",{
  #select samples that are in all three types of mutations
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                                          gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- unique(samples)[1:5]

  sp <- as.data.frame(samples) %>%
    mutate(panels = "IMPACT341")

  mut <- gnomeR::mutations %>%
    filter(sampleId %in% samples)


  sp <- rename(sp,
               "sample_id" = samples,
               "panel_id" = panels)

  sp[1, 2] <- NA_character_

  expect_error(bin_impact <-  create_gene_binary(samples=samples,
                                                   mutation = mut,
                                                   cna = gnomeR::cna,
                                                   fusion = gnomeR::sv,
                                                   specify_panel = sp), "*")




})


