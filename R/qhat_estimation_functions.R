#' Estimate marginal and joint distribution of lists i and j using logistic regression.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_logit(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_logit <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, ...){

  fiti0 = try(glm(formula(paste("L", i, "*(1 - L", j, ") ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(i, j, (K + 1):ncol(List.train))]))
  fit0j = try(glm(formula(paste("L", j, "*(1 - L", i, ") ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(i, j, (K + 1):ncol(List.train))]))
  fitij = try(glm(formula(paste("L", i, "*L", j, " ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(i, j, (K + 1):ncol(List.train))]))

  if("try-error" %in% c(class(fiti0), class(fit0j), class(fitij))){
    Warning("One or more fits with logistic regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitij, newdata = List.test, type = "response"), eps)
    q1 = pmax(q12 + predict(fiti0, newdata = List.test, type = "response"), eps)
    q2 = pmax(q12 + predict(fit0j, newdata = List.test, type = "response"), eps)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}

#' Estimate marginal and joint distribution of lists i and j using generalized additive models.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_gam(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_gam <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, ...){

  require("gam")
  l = ncol(List.train) - K
  #colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))
  #colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))

  #sapply((K + 1):ncol(List.train), function(l){ paste(rep(paste('x', l, sep = ' + '), 3), c(' + ', '^2 + ', '^3 + '), sep = '') })
  fiti = try(gam::gam(formula(paste("L", i, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List.train, family = binomial("logit")))
  fitj = try(gam::gam(formula(paste("L", j, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List.train, family = binomial("logit")))
  fitij = try(gam::gam(formula(paste("L", i, "*L", j, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List.train, family = binomial("logit")))

  if("try-error" %in% c(class(fiti), class(fitj), class(fitij))){
    Warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitij, newdata = List.test, type = "response"), eps)
    q1 = pmax(predict(fiti, newdata = List.test, type = "response"), eps)
    q2 = pmax(predict(fitj, newdata = List.test, type = "response"), eps)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}

#' Estimate marginal and joint distribution of lists i and j using super learner.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_sl(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_sl <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet")){

  require("SuperLearner")

  slib = intersect(sl.lib, c("SL.glm", "SL.gam", "SL.glm.interaction"))
  slib1 = setdiff(sl.lib, slib)

  slib2 <- c(slib1, slib,
             split(rbind(slib,"screen.corP"),
                   rep(1:length(slib),each=2)),
             split(rbind(slib,"screen.glmnet"),
                   rep(1:length(slib),each=2)))

  if(ncol(List.train) == K + 1){
    xtrain = data.frame(x1 = List.train[,K+1])
    xtest = data.frame(x1 = List.test[,K+1])
  }else{
    xtrain = List.train[,-c(1:K)]
    xtest = List.test[,-c(1:K)]
  }

  fiti = tryCatch(SuperLearner(Y = as.numeric(List.train[,i]),
                          X = xtrain,
                          family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)
  fitj = tryCatch(SuperLearner(Y = as.numeric(List.train[,j]),
                          X = xtrain,
                          family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)
  fitij = tryCatch(SuperLearner(Y = as.numeric(pmin(List.train[,i], List.train[,j])),
                           X = xtrain,
                           family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)

  if("try-error" %in% c(class(fiti), class(fitj), class(fitij))){
    Warning("One or more fits with SuperLearner regression failed.")
    return(NULL)
  }else{
    q12 = pmax(pmin(
      predict(fitij, newdata = xtest, onlySL = TRUE)$pred, 1), eps)
    q1 = pmin(pmax(
      predict(fiti, newdata = xtest, onlySL = TRUE)$pred, q12), 1)
    q2 = pmin(pmax(
      predict(fitj, newdata = xtest, onlySL = TRUE)$pred, q12), 1)
  }

  return(list(q1 = q1, q2 = q2, q12 = q12))
}
