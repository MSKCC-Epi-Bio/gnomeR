
#setting these patients as global variables. so that we can use them in all tests throughout
set.seed(123)
patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
