# context("check make maf summary")
# source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())
#
# test_that("binary cna",{
#   expect_warning(gen_dat <- binmat(samples = samples, maf = mut, cna = cna, cna.binary = TRUE))
#   expect_warning(test <- gen_uni_cox(X = gen_dat,surv_dat = surv_dat,surv_formula  = Surv(time,status)~.,filter = 0,genes = NULL))
#   expect_true(is.data.frame(test$tab))
#   expect_true(nrow(test$tab) == 613)
#   expect_true(all(class(test$p) == c("plotly","htmlwidget")))
#   expect_true(all(class(test$KM[[1]]) == c("ggsurvplot","ggsurv","list")))
# })
#
# test_that("non binary cna",{
#   expect_warning(gen_dat <- binmat(samples = samples, maf = mut, cna = cna, cna.binary = FALSE))
#   expect_warning(test <- gen_uni_cox(X = gen_dat,surv_dat = surv_dat,surv_formula  = Surv(time,status)~.,filter = 0,genes = NULL))
#   expect_true(nrow(test$tab) == 624)
#   expect_true(all(class(test$p) == c("plotly","htmlwidget")))
#   expect_true(all(class(test$KM[[1]]) == c("ggsurvplot","ggsurv","list")))
# })
