library(cbioportalR)
library(dplyr)

set_cbioportal_db("public")
#save <- available_studies()
gen <- get_genetics_by_study("prad_mskcc_2017")


mutations <- gen$mutation
cna <- gen$cna
sv <- gen$structural_variant

#sample 200 patients
un <-  unique(mut$sampleId)
sample_patients <- sample(un, size = 200, replace = FALSE)

#subset the three datasets for those sampled patients
mutations <- mutations %>%
  filter(sampleId %in% sample_patients)

cna <- cna %>%
  filter(sampleId %in% sample_patients)

sv <- sv %>%
  filter(sampleId %in% sample_patients)



usethis::use_data(mutations, overwrite = TRUE)
usethis::use_data(cna, overwrite = TRUE)
usethis::use_data(sv, overwrite = TRUE)
