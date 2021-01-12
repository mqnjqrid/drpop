#' Estimate marginal and joint distribution of lists i and j using logistic regression.
#'
#' @param List1 The training data matrix used to estimate the distibution functions.
#' @param List2 The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_logit(List1 = List1, List2 = List2, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_logit <- function(List1, List2, K, i, j, eps){

  fiti0 = try(glm(formula(paste("L", i, "*(1 - L", j, ") ~.", sep = '')), family = binomial(link = "logit"), data = List1[,c(i, j, (K + 1):ncol(List1))]))
  fit0j = try(glm(formula(paste("L", j, "*(1 - L", i, ") ~.", sep = '')), family = binomial(link = "logit"), data = List1[,c(i, j, (K + 1):ncol(List1))]))
  fitij = try(glm(formula(paste("L", i, "*L", j, " ~.", sep = '')), family = binomial(link = "logit"), data = List1[,c(i, j, (K + 1):ncol(List1))]))

  if("try_error" %in% c(class(fiti0), class(fit0j), class(fitij))){
    Warning("One or more fits with logistic regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitij, newdata = List2, type = "response"), eps)
    q1 = pmax(q12 + predict(fiti0, newdata = List2, type = "response"), eps)
    q2 = pmax(q12 + predict(fit0j, newdata = List2, type = "response"), eps)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}

#' Estimate marginal and joint distribution of lists i and j using generalized additive models.
#'
#' @param List1 The training data matrix used to estimate the distibution functions.
#' @param List2 The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_gam(List1 = List1, List2 = List2, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_gam = function(List1, List2, K, i, j, eps){

  l = ncol(List1) - K
  #colnames(List1) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))
  #colnames(List2) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))

  #sapply((K + 1):ncol(List1), function(l){ paste(rep(paste('x', l, sep = ' + '), 3), c(' + ', '^2 + ', '^3 + '), sep = '') })
  fiti = try(gam::gam(formula(paste("L", i, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List1, family = binomial))
  fitj = try(gam::gam(formula(paste("L", j, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List1, family = binomial))
  fitij = try(gam::gam(formula(paste("L", i, "*L", j, " ~", paste("s(x", 1:l, ")", sep = '', collapse = ' + '), sep = '')), data = List1, family = binomial))

  if("try_error" %in% c(class(fiti0), class(fit0j), class(fitij))){
    Warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitij, newdata = List2, type = "response"), eps)
    q1 = pmax(q12 + predict(fiti, newdata = List2, type = "response"), eps)
    q2 = pmax(q12 + predict(fitj, newdata = List2, type = "response"), eps)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}

#' Estimate marginal and joint distribution of lists i and j using super learner.
#'
#' @param List1 The training data matrix used to estimate the distibution functions.
#' @param List2 The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_sl(List1 = List1, List2 = List2, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_sl = function(List1, List2, K, i, j, eps){
  slib = c("SL.glm"
           , "SL.gam2"
           , "SL.glm.interaction"
  )
  slib1 = c("SL.glmnet"
            ,"SL.ranger"
            #, "SL.gbm"
  )
  slib2 <- c(slib1, slib,
             split(rbind(slib,"screen.corP"),
                   rep(1:length(slib),each=2)) ,
             split(rbind(slib,"screen.glmnet"),
                   rep(1:length(slib),each=2)) )

  #suppressWarnings()

  fiti = try(SuperLearner(Y = as.numeric(List1[,i]),
                          X = as.data.frame(List1[,-c(1:K)]),
                          family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)
  fitj = try(SuperLearner(Y = as.numeric(List1[,j]),
                          X = as.data.frame(List1[,-c(1:K)]),
                          family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)
  fitij = try(SuperLearner(Y = as.numeric(pmin(List1[,i], List1[,j])),
                           X = as.data.frame(List1[,-c(1:K)]),
                           family = binomial(), SL.library = slib2, verbose = FALSE), silent = TRUE)

  if("try_error" %in% c(class(fiti0), class(fit0j), class(fitij))){
    Warning("One or more fits with SuperLearner regression failed.")
    return(NULL)
  }else{
    q12 = pmax(pmin(
      predict(fitij, newdata = List2[,-c(1:K)], onlySL = TRUE)$pred, 1), eps)
    q1 = pmin(pmax(
      predict(fiti, newdata = List2[,-c(1:K)], onlySL = TRUE)$pred, q12), 1)
    q2 = pmin(pmax(
      predict(fitj, newdata = List2[,-c(1:K)], onlySL = TRUE)$pred, q12), 1)
  }

  return(list(q1 = q1, q2 = q2, q12 = q12))
}
