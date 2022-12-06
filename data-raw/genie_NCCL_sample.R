library(genieBPC)
library(dplyr)


nsclc_2_0 = pull_data_synapse("NSCLC", version = "v2.0-public")


genie_mut <- nsclc_2_0$NSCLC_v2.0[[11]]
genie_cna <- nsclc_2_0$NSCLC_v2.0[[10]]
genie_fusion <- nsclc_2_0$NSCLC_v2.0[[9]]

#sample 200 patients
un <-  unique(genie_mut$Tumor_Sample_Barcode)
sample_patients <- sample(un, size = 100, replace = FALSE)

#for cna need to subset the columns not rows
cna_samps <- gsub("-", ".",sample_patients)
cna_samps2 <- colnames(genie_cna)[colnames(genie_cna) %in% cna_samps]

#subset the three datasets for those sampled patients
genie_mut <- genie_mut %>%
  filter(Tumor_Sample_Barcode %in% sample_patients)

genie_cna <- genie_cna %>%
  select(Hugo_Symbol, all_of(cna_samps2))

genie_fusion <- genie_fusion %>%
  filter(Tumor_Sample_Barcode %in% sample_patients)



usethis::use_data(genie_mut, overwrite = TRUE)
usethis::use_data(genie_cna, overwrite = TRUE)
usethis::use_data(genie_fusion, overwrite = TRUE)
