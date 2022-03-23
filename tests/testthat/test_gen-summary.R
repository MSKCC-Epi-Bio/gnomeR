context("check make gen summary")
source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())

test_that("wrong filter",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "A","B")
  outcome <- sample(c("A","B"),nrow(mat),replace = T)
  expect_error(gen.summary(gen.dat = mat,
                       outcome = outcome,
                       filter = -1))
  expect_error(gen.summary(gen.dat = mat,
                       outcome = outcome,
                       filter = 1))
  expect_error(gen.summary(gen.dat = mat,
                       outcome = outcome,
                       filter = 2))

})


test_that("working binary example",{

  outcome <- sample(c("A","B"),nrow(gen.dat),replace = T)
  test <- gen.summary(gen.dat = gen.dat,
  outcome = outcome,
  filter = 0.05,
  cont = FALSE,rank = TRUE)

  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 9)
  expect_true(is.ggplot(test$forest.plot))
  # expect_true(is.ggplot(test$vPlot))
})


test_that("working continuous example",{

  gen.dat.test <- gen.dat
  gen.dat.test[,1] <- rnorm(n = nrow(gen.dat.test))
  outcome <-  rnorm(n = nrow(gen.dat.test))
  test <- gen.summary(gen.dat = gen.dat.test,
                     outcome = outcome,
                     filter = 0.05,
                     cont = TRUE,rank = TRUE)


  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 8)
  expect_true(is.null(test$forest.plot))
  # expect_true(!is.null(test$vPlot))
})


test_that("working binary example with a feature that is only 1's",{

  gen.dat.test<- gen.dat
  gen.dat.test[,1] <- 1
  outcome <-  rnorm(n = nrow(gen.dat.test))
  test <- gen.summary(gen.dat = gen.dat.test,
                  outcome = outcome,
                  filter = 0.05,
                  cont = TRUE,rank = TRUE)

  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 8)
  expect_true(all(class(test$vPlot) == c("plotly","htmlwidget")))
})


test_that("filter too large",{

  gen.dat.test<-gen.dat
  expect_error(gen.summary(gen.dat = gen.dat.test,
                  outcome = outcome,
                  filter = 0.99,
                  cont = FALSE,rank = TRUE))

})

# test_that("continuous features",{
#
#   gen.dat.test<-gen.dat
#   gen.dat.test[,1] <- rnorm(n = nrow(gen.dat.test))
#   test <- gen.summary(gen.dat = gen.dat.test,
#           outcome = outcome,
#           filter = 0,
#           cont = FALSE,rank = TRUE)
#   expect_true(is.data.frame(test$fits))
#   expect_true(ncol(test$fits) == 9)
#   expect_true(is.ggplot(test$forest.plot))
#   # expect_true(is.ggplot(test$vPlot))
#
# })


# test_that("paired test",{
#
#   set.seed(123)
#   gen.dat <- as.data.frame(matrix(rbinom(500,1,1/2),nrow = 100, ncol = 5))
#   outcome <- c(rep("Time1",50),rep("Time2",50))
#   test <- gen.summary(gen.dat = gen.dat,
#                   outcome = outcome,
#                   filter = 0,paired = TRUE,
#                   cont = FALSE,rank = TRUE)
#   expect_true(is.data.frame(test$fits))
#   expect_true(ncol(test$fits) == 9)
#   expect_true(is.ggplot(test$forest.plot))
#   # expect_true(is.ggplot(test$vPlot))
#
# })


test_that("three level outcome",{

  set.seed(123)
  gen.dat <- as.data.frame(matrix(rbinom(500,1,1/2),nrow = 100, ncol = 5))
  outcome <- sample(c("A","B","C"),100,replace = TRUE)
  test <- gen.summary(gen.dat = gen.dat,
                  outcome = outcome,
                  filter = 0,
                  cont = FALSE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 10)
  expect_true(is.null(test$forest.plot))
  expect_true(is.null(test$vPlot))

})

test_that("factors hidden in continuous variables",{

  gen.dat.test <- gen.dat
  gen.dat.test[,1] <- factor(sample(c("NEUTRAL","DELETION","LOH"),nrow(gen.dat.test),replace = TRUE),
                             levels = c("NEUTRAL","LOH","DELETION"))
  outcome <-  rnorm(n = nrow(gen.dat.test))
  test <- gen.summary(gen.dat = gen.dat.test[,-1],
                  outcome = outcome,
                  filter = 0,
                  cont = TRUE,rank = TRUE)
  expect_true(is.data.frame(test$fits))
  expect_true(ncol(test$fits) == 8)
  expect_true(is.null(test$forest.plot))
  expect_true(!is.null(test$vPlot))

})

