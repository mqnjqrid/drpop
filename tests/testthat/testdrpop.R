test_that("qhat_logit estimation runs fine", {
  data = simuldata(n = 2000, l = 1)$data
  qlogit = qhat_logit(List.train = data)
  expect_false(class(qlogit) == "try-error")
})
test_that("qhat_gam estimation runs fine", {
  data = simuldata(n = 1000, l = 1)$data
  qgam = qhat_gam(List.train = data)
  expect_false(class(qgam) == "try-error")
})

test_that("qhat_sl estimation runs fine", {
  data = simuldata(n = 1000, l = 3)$data
  qsl = qhat_sl(List.train = data)
  expect_false(class(qsl) == "try-error")
})

test_that("qhat_ranger estimation runs fine", {
  data = simuldata(n = 1000, l = 1)$data
  qranger = qhat_ranger(List.train = data)
  expect_false(class(qranger) == "try-error")
})

test_that("Runs for single covariate", {
  ps = try(drpop::popsize(data = simuldata(n = 1000, l = 1)$data, funcname = c("logit", "gam", "mlogit")), silent = TRUE)
  expect_false(class(ps) == "try-error")
  expect_equal(sum(sapply(ps[-6], anyNA)), 0)
})

test_that("Runs for single covariate with factor column", {
  ps = try(popsize(data = simuldata(n = 1000, l = 1, categorical = TRUE)$data, K = 2, funcname = c("logit", "gam", "mlogit", "ranger")), silent = TRUE)
  expect_false(class(ps) == "try-error")
  expect_equal(sum(sapply(ps[-6], anyNA)), 0)
})

#sdata = simuldata(n = 1000, l = 3, categorical = FALSE)
#data = sdata$data

#ps = popsize(data = data)

sdata = simuldata(n = 1000, l = 3, categorical = TRUE)
data = sdata$data
ps = popsize(data = data)
psq = popsize(data = data, getnuis = ps$nuis, idfold = ps$idfold)
test_that("Returns no NA values", {
  expect_equal(sum(sapply(ps$result[-6], anyNA)), 0)
  expect_equal(sum(sapply(psq$result[-6], anyNA)), 0)
})

test_that("Error catching",{
  expect_error(popsize(data = c(1:10)))
  expect_error(popsize(data = simuldata(n = 7000, l = 6)$data, K = 5))
})

test_that("Subset of lists", {
  ps = try(popsize(data = simuldata(n = 1000, l = 3, K = 3)$data, K = 3, j = 2, k = 3, funcname = c("gam", "logit")))
  expect_equal(as.character(unique(ps$result$listpair)), c("2,3"))
})

test_that("All lists", {
  ps = try(popsize(data = simuldata(n = 1000, l = 3, K = 3)$data, K = 3, funcname = c("gam", "logit")))
  expect_setequal(as.character(unique(ps$result$listpair)), c( "1,2", "2,3", "1,3"))
})
