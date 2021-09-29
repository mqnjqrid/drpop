
test_that("Runs for single covariate", {
  ps = try(drpop::popsize(List_matrix = simuldata(n = 1000, l = 1)$List_matrix, funcname = c("logit", "gam", "mlogit")), silent = TRUE)
  print(ps)
  expect_equal(class(ps) == "try-error", FALSE)
  expect_equal(sum(sapply(ps, anyNA)), 0)
})

test_that("Runs for single covariate with factor column", {
  ps = try(popsize(List_matrix = simuldata(n = 1000, l = 1, categorical = TRUE)$List_matrix, funcname = c("logit", "gam", "mlogit", "ranger")), silent = TRUE)
  expect_equal(class(ps) == "try-error", FALSE)
  expect_equal(sum(sapply(ps, anyNA)), 0)
})
sdata = simuldata(n = 1000, l = 3, categorical = FALSE)
List_matrix = sdata$List_matrix

ps = popsize(List_matrix = List_matrix)

sdata = simuldata(n = 1000, l = 3, categorical = TRUE)
List_matrix = sdata$List_matrix
ps = popsize(List_matrix = List_matrix)
psq = popsize(List_matrix = List_matrix, getnuis = ps$nuis, idfold = ps$idfold)
test_that("Returns no NA values", {
  expect_equal(sum(sapply(ps$result, anyNA)), 0)
  expect_equal(sum(sapply(psq$result, anyNA)), 0)
})

test_that("Error catching",{
  expect_error(popsize(List_matrix = c(1:10)))
  expect_error(popsize(List_matrix = simuldata(n = 7000, l = 6)$List_matrix, K = 5))
})

test_that("Subset of lists", {
  ps = try(popsize(List_matrix = simuldata(n = 1000, l = 3, K = 3)$List_matrix, K = 3, j = 2, k = 3, funcname = c("gam", "logit")))
  expect_equal(as.character(unique(ps$result$listpair)), c("2,3"))
})

test_that("All lists", {
  ps = try(popsize(List_matrix = simuldata(n = 1000, l = 3, K = 3)$List_matrix, K = 3, funcname = c("gam", "logit")))
  expect_setequal(as.character(unique(ps$result$listpair)), c( "1,2", "2,3", "1,3"))
})
