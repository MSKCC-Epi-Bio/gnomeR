
#setting these patients as global variables. so that we can use them in all tests throughout
set.seed(123)
patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 300, replace=FALSE)]
time <- rexp(300)
status <- rbinom(n = 300,size = 1,prob = 1/2)
surv.dat <- as.data.frame(cbind(time,status))

#global variables for gen-sumamry
patients100 <- patients[1:100]
outcome <- as.character(clin.sample$Sample.Type[match(patients100,clin.sample$Sample.Identifier)])
gen.dat <- try(binmat(patients = patients100,maf = mut),silent = TRUE)
