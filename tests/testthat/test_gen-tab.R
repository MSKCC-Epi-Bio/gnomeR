context("check make gen tab")

test_that("wrong filter",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "A","B")
  outcome <- sample(c("A","B"),nrow(mat),replace = T)
  expect_error(gen.tab(gen.dat = mat,
                       outcome = outcome,
                       filter = -1))
  expect_error(gen.tab(gen.dat = mat,
                       outcome = outcome,
                       filter = 1))
  expect_error(gen.tab(gen.dat = mat,
                       outcome = outcome,
                       filter = 2))

})


test_that("working binary example",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  outcome <- as.character(clin.sample$Sample.Type[match(patients,clin.sample$Sample.Identifier)])
  gen.dat <- binmat(patients = patients,maf = mut)
  test <- gen.tab(gen.dat = gen.dat,
  outcome = outcome,
  filter = 0.05,paired = FALSE,
  cont = FALSE,rank = TRUE)

  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 9)
  expect_true(is.ggplot(test$forest.plot))
  expect_true(is.ggplot(test$vPlot))
})


test_that("working continuous example",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  gen.dat <- binmat(patients = patients,maf = mut)
  gen.dat[,1] <- rnorm(n = nrow(gen.dat))
  outcome <-  rnorm(n = nrow(gen.dat))
  test <- gen.tab(gen.dat = gen.dat,
                     outcome = outcome,
                     filter = 0.05,paired = FALSE,
                     cont = TRUE,rank = TRUE)


  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 5)
  expect_true(is.null(test$forest.plot))
  expect_true(!is.null(test$vPlot))
})


test_that("working binary example with a feature that is only 1's",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  outcome <- as.character(clin.sample$Sample.Type[match(patients,clin.sample$Sample.Identifier)])
  gen.dat <- binmat(patients = patients,maf = mut)
  gen.dat[,1] <- 1
  test <- gen.tab(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0.05,paired = FALSE,
                  cont = FALSE,rank = TRUE)

  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 9)
  expect_true(is.ggplot(test$forest.plot))
  expect_true(is.ggplot(test$vPlot))
})


test_that("filter too large",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  outcome <- as.character(clin.sample$Sample.Type[match(patients,clin.sample$Sample.Identifier)])
  gen.dat <- binmat(patients = patients,maf = mut)
  expect_error(gen.tab(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0.99,paired = FALSE,
                  cont = FALSE,rank = TRUE))

})

test_that("continuous features",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  outcome <- as.character(clin.sample$Sample.Type[match(patients,clin.sample$Sample.Identifier)])
  gen.dat <- binmat(patients = patients,maf = mut)
  gen.dat[,1] <- rnorm(n = nrow(gen.dat))
  test <- gen.tab(gen.dat = gen.dat,
          outcome = outcome,
          filter = 0,paired = FALSE,
          cont = FALSE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 9)
  expect_true(is.ggplot(test$forest.plot))
  expect_true(is.ggplot(test$vPlot))

})


test_that("paired test",{

  set.seed(123)
  gen.dat <- as.data.frame(matrix(rbinom(500,1,1/2),nrow = 100, ncol = 5))
  outcome <- c(rep("Time1",50),rep("Time2",50))
  test <- gen.tab(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0,paired = TRUE,
                  cont = FALSE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 9)
  expect_true(is.ggplot(test$forest.plot))
  expect_true(is.ggplot(test$vPlot))

})


test_that("three level outcome",{

  set.seed(123)
  gen.dat <- as.data.frame(matrix(rbinom(500,1,1/2),nrow = 100, ncol = 5))
  outcome <- sample(c("A","B","C"),100,replace = TRUE)
  test <- gen.tab(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0,paired = FALSE,
                  cont = FALSE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 10)
  expect_true(is.null(test$forest.plot))
  expect_true(is.null(test$vPlot))

})

test_that("factors hidden in continuous variables",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:100]
  gen.dat <- binmat(patients = patients,maf = mut)
  gen.dat[,1] <- sample(c("0","-2","-2"),nrow(gen.dat),replace = TRUE)
  outcome <-  rnorm(n = nrow(gen.dat))
  test <- gen.tab(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0,paired = FALSE,
                  cont = TRUE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 5)
  expect_true(is.null(test$forest.plot))
  expect_true(!is.null(test$vPlot))

})

