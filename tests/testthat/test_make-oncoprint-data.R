context("check make oncoprint dat")

test_that("make data nothing special",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut,mut.type = "SOMATIC",
  SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
  gen.dat <- bin.mut[,
  names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
  out <- dat.oncoPrint(gen.dat)
  expect_true(is.matrix(out))
  expect_true(all(dim(out) == c(15,200)))

})


test_that("make data while throwing warnings",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut,mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
  gen.dat <- bin.mut[,
                     names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
  gen.dat[,1] <- sample(c("A","B"), size = nrow(gen.dat), replace = TRUE)
  expect_warning(out <- dat.oncoPrint(gen.dat))
  expect_true(is.matrix(out))
  expect_true(all(dim(out) == c(14,200)))

})


test_that("make data while throwing warnings",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
  bin.mut <- binmat(patients = patients,maf = mut,mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
  gen.dat <- bin.mut[,
                     names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
  gen.dat[,1] <- sample(c(0:2,NA), size = nrow(gen.dat), replace = TRUE)
  expect_warning(out <- dat.oncoPrint(gen.dat))
  expect_true(is.matrix(out))
  expect_true(all(dim(out) == c(14,200)))

})


test_that("make data with clinical",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))
  bin.mut <- binmat(patients = patients,maf = mut, fusion = fusion,
                    cna = cna, mut.type = "SOMATIC",
                    SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
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
    dplyr::mutate(Smoker = ifelse(Smoker == "Unknown",NA, as.character(Smoker)))
  clin.patients.dat$ContFeature <- rnorm(n = nrow(clin.patients.dat))
  gen.dat <- gen.dat[match(clin.patients.dat$DMPID,
  abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9)),]
  clin.patients.dat <- clin.patients.dat %>%
    tibble::column_to_rownames('DMPID')
  rownames(gen.dat) <- rownames(clin.patients.dat)
  out <- dat.oncoPrint(gen.dat,clin.patients.dat)
  expect_true(is.matrix(out))
  expect_true(all(dim(out) == c(6,454)))

})
