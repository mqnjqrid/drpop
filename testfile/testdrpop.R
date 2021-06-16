library(testthat)

devtools::load_all("C:/Users/manja/OneDrive/Documents/drpop")
test_that("Runs for single covariate", {
  ps = try(popsize(List_matrix = simuldata(n = 1000, l = 1)$List_matrix, funcname = c("rangerlogit", "gam", "mlogit", "ranger")))
  expect_equal(class(ps) == "try-error", FALSE)
  expect_equal(sum(sapply(ps, anyNA)), 0)
})

test_that("Runs for single covariate with factor column", {
  ps = try(popsize(List_matrix = simuldata(n = 1000, l = 1, categorical = TRUE)$List_matrix, funcname = c("logit", "gam", "mlogit", "ranger")))
  expect_equal(class(ps) == "try-error", FALSE)
  expect_equal(sum(sapply(ps, anyNA)), 0)
})
sdata = simuldata(n = 1000, l = 3, categorical = FALSE)
List_matrix = sdata$List_matrix

ps = popsize(List_matrix = List_matrix)

qhatev = getnuis(List_matrix = List_matrix)
psq = popsizenuis(List_matrix = List_matrix, getnuis = qhatev)

test_that("Returns no NA values", {
  expect_equal(sum(sapply(ps, anyNA)), 0)
  expect_equal(sum(sapply(psq, anyNA)), 0)
})

sdata = simuldata(n = 1000, l = 3, categorical = TRUE)
List_matrix = sdata$List_matrix
ps = popsize(List_matrix = List_matrix)
qhatev = getnuis(List_matrix = List_matrix)
psq = popsizenuis(List_matrix = List_matrix, getnuis = qhatev)
test_that("Returns no NA values", {
  expect_equal(sum(sapply(ps, anyNA)), 0)
  expect_equal(sum(sapply(psq, anyNA)), 0)
})

test_that("Error catching",{
  expect_error(popsizenuis(List_matrix = c(1:10)))
  expect_error(popsizenuis(List_matrix = simuldata(n = 7000, l = 6)$List_matrix))
})
