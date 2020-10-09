context("check make oncoprint")


test_that("make data genetics only",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut, fusion = fusion,
                    cna = cna, mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, specify.plat = FALSE)
  gen.dat <- bin.mut[,
                     c("TP53","CDKN2A","CDKN2A.Del","PIK3CA","PIK3CA.Amp")]
  gen.dat$TP53 <- rbinom(nrow(gen.dat),1,1/2)
  gen.dat$CDKN2A.Del <- rbinom(nrow(gen.dat),1,1/2)
  out <- plot_oncoPrint(gen.dat)
  expect_true(typeof(out) == "S4")
  out <- plot_oncoPrint(gen.dat,ordered = 1:nrow(gen.dat))
  expect_true(typeof(out) == "S4")

})

test_that("make data with clinical",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut, fusion = fusion,
                    cna = cna, mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, specify.plat = FALSE)
  gen.dat <- bin.mut[,
                     c("TP53", "CDKN2A","CDKN2A.Del","PIK3CA","PIK3CA.Amp")]
  gen.dat$TP53 <- rbinom(nrow(gen.dat),1,1/2)
  gen.dat$CDKN2A.Del <- rbinom(nrow(gen.dat),1,1/2)
  clin.patients.dat <-
    clin.patients[match(abbreviate(rownames(gen.dat),
                                   strict = TRUE, minlength = 9),
                        clin.patients$X.Patient.Identifier),] %>%
    dplyr::rename(DMPID = X.Patient.Identifier,
                  Smoker = Smoking.History) %>%
    dplyr::select(DMPID, Sex,Smoker) %>%
    dplyr::filter(!is.na(DMPID)) %>%
    dplyr::distinct(DMPID,.keep_all = TRUE) %>%
    dplyr::mutate(Smoker = ifelse(Smoker == "Unknown",NA, as.character(Smoker)))
  clin.patients.dat$ContFeature <- rnorm(n = nrow(clin.patients.dat))
  gen.dat <- gen.dat[match(clin.patients.dat$DMPID,
                           abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9)),]
  clin.patients.dat <- clin.patients.dat %>%
    tibble::column_to_rownames('DMPID')
  rownames(gen.dat) <- rownames(clin.patients.dat)
  out <- plot_oncoPrint(gen.dat,clin.patients.dat,ordered = 1:nrow(gen.dat))
  expect_true(typeof(out) == "S4")

})


test_that("make data single clinical",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut, fusion = fusion,
                    cna = cna, mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, specify.plat = FALSE)
  gen.dat <- bin.mut[,
                     c("TP53","CDKN2A","CDKN2A.Del","PIK3CA","PIK3CA.Amp")]
  gen.dat$TP53 <- rbinom(nrow(gen.dat),1,1/2)
  gen.dat$CDKN2A.Del <- rbinom(nrow(gen.dat),1,1/2)
  clin.patients.dat <-
    clin.patients[match(abbreviate(rownames(gen.dat),
                                   strict = TRUE, minlength = 9),
                        clin.patients$X.Patient.Identifier),] %>%
    dplyr::rename(DMPID = X.Patient.Identifier,
                  Smoker = Smoking.History) %>%
    dplyr::select(DMPID, Sex,Smoker) %>%
    dplyr::filter(!is.na(DMPID)) %>%
    dplyr::distinct(DMPID,.keep_all = TRUE) %>%
    dplyr::mutate(Smoker = ifelse(Smoker == "Unknown",NA, as.character(Smoker))) %>%
    dplyr::select(Sex,DMPID)
  clin.patients.dat$ContFeature <- rnorm(n = nrow(clin.patients.dat))
  clin.patients.dat <- clin.patients.dat %>%
    dplyr::select(ContFeature,DMPID)
  gen.dat <- gen.dat[match(clin.patients.dat$DMPID,
                           abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9)),]
  clin.patients.dat <- clin.patients.dat %>%
    tibble::column_to_rownames('DMPID')
  rownames(gen.dat) <- rownames(clin.patients.dat)
  out <- plot_oncoPrint(gen.dat,clin.dat = clin.patients.dat)
  expect_true(typeof(out) == "S4")

})


test_that("make data with cna non binary",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut,cna = cna, cna.binary = FALSE,
                    mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, specify.plat = FALSE)
  keep <- c("TP53|PIK3CA|ALK")
  gen.dat <- bin.mut[,grep(keep,colnames(bin.mut))]
  expect_warning(out <- plot_oncoPrint(gen.dat))
  expect_true(typeof(out) == "S4")

})

