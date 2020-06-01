context("check make maf summary")
source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())

test_that("binary cna",{
  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 300, replace=FALSE)]
  gen.dat <- binmat(patients = patients, maf = mut, cna = cna, cna.binary = TRUE)
  time <- rexp(nrow(gen.dat))
  status <- rbinom(n = nrow(gen.dat),size = 1,prob = 1/2)
  surv.dat <- as.data.frame(cbind(time,status))
  expect_warning(test <- uni.cox(X = gen.dat,surv.dat = surv.dat,surv.formula  = Surv(time,status)~.,filter = 0,genes = NULL))
  expect_true(is.data.frame(test$tab))
  expect_true(nrow(test$tab) == 603)
  expect_true(all(class(test$p) == c("plotly","htmlwidget")))
  expect_true(all(class(test$KM[[1]]) == c("ggsurvplot","ggsurv","list")))
})

test_that("non binary cna",{
  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 300, replace=FALSE)]
  gen.dat <- binmat(patients = patients, maf = mut, cna = cna, cna.binary = FALSE)
  time <- rexp(nrow(gen.dat))
  status <- rbinom(n = nrow(gen.dat),size = 1,prob = 1/2)
  surv.dat <- as.data.frame(cbind(time,status))
  test <- uni.cox(X = gen.dat,surv.dat = surv.dat,surv.formula  = Surv(time,status)~.,filter = 0,genes = NULL)
  expect_true(nrow(test$tab) == 644)
  expect_true(all(class(test$p) == c("plotly","htmlwidget")))
  expect_true(all(class(test$KM[[1]]) == c("ggsurvplot","ggsurv","list")))
})
