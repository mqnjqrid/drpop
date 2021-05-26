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

  stopifnot(ncol(List.train) > K)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

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

  stopifnot(ncol(List.train) > K)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  require("gam", quietly = TRUE)
  l = ncol(List.train) - K
  nnum = intersect(names(which(lapply(subset(List.train, select = -c(1:K)), "class") != "factor")),
                 names(which(lapply(subset(List.train, select =  -c(1:K)), function(icol){length(unique(icol))}) > 3)))
  nfac = setdiff(names(List.train[,-c(1:K)]), nnum)
  if(!length(nnum)){
    covformula = paste0(nfac, collapse = ' + ')
  }else if(!length(nfac)){
    covformula = paste0('s(', nnum, ')', collapse = ' + ')
  }else{
    covformula = paste0(paste0('s(', nnum, ')', collapse = ' + '),
                        '+', paste0('s(', nnum, ')', collapse = ' + '))
  }

  fiti = try(gam::gam(formula(paste0("L", i, " ~", covformula)), data = List.train, family = binomial("logit")))
  fitj = try(gam::gam(formula(paste0("L", j, " ~", covformula)), data = List.train, family = binomial("logit")))
  fitij = try(gam::gam(formula(paste0("L", i, "*L", j, " ~", covformula)), data = List.train, family = binomial("logit")))

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
#' Estimate marginal and joint distribution of lists i and j using ranger.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_ranger(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_ranger <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, ...){
  require("ranger", quietly = TRUE)
  l = ncol(List.train) - K

  stopifnot(l>0)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  fiti = ranger(formula(paste("factor(L", i, ") ~.", sep = '')), data = List.train[,-c(1:K)[-i]], probability = TRUE, classification = TRUE)
  fitj = ranger(formula(paste("factor(L", j, ") ~.", sep = '')), data = List.train[,-c(1:K)[-j]], probability = TRUE, classification = TRUE)
  fitij = ranger(formula(paste("factor(L", i, "*L", j, ") ~.", sep = '')), data = List.train[,c(i, j, K + 1:l)], probability = TRUE, classification = TRUE)

  if("try-error" %in% c(class(fiti), class(fitj), class(fitij))){
    Warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitij, data = List.test, type = "response")$predictions[,'1'], eps)
    q1 = pmax(predict(fiti, data = List.test, type = "response")$predictions[,'1'], eps)
    q2 = pmax(predict(fitj, data = List.test, type = "response")$predictions[,'1'], eps)

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
#' @param sl.lib The functions from the SuperLearner library to be used for model fitting.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param num_cores The number of cores to be used for paralellization in Super Learner.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_sl(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_sl <- function (List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005,
                     sl.lib = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), num_cores = NA,...)
{

  stopifnot(ncol(List.train) > K)

  require("SuperLearner", quietly = TRUE)
  require("parallel", quietly = TRUE)
  require("gam", quietly = TRUE)
  require("xgboost", quietly = TRUE)
  require("janitor", quietly = TRUE)
  require("tidyverse", quietly = TRUE)
  slib = intersect(sl.lib, c("SL.glm", "SL.gam",
                             "SL.glm.interaction"))
  slib1 = setdiff(sl.lib, slib)
  if(length(slib) == 0){
    slib2 = slib1
  }else{
    slib2 <- c(slib1, slib, split(rbind(slib, "screen.corP"),
                                  rep(1:length(slib), each = 2)), split(rbind(slib, "screen.glmnet"),
                                                                        rep(1:length(slib), each = 2)))
  }
  num_cores = min(num_cores, parallel::detectCores() - 1, na.rm = TRUE)

  factor_cols <- subset(List.train, select = -c(1:K)) %>% select_if(negate(is.numeric)) %>% names()

  if(length(factor_cols)) {
    List.train = data.frame(List.train[,!(names(List.train) %in% factor_cols)],
                              model.matrix(formula(paste('~', paste(factor_cols, collapse = '+'))), List.train))
    List.test = data.frame(List.test[,!(names(List.test) %in% factor_cols)],
                            model.matrix(formula(paste('~', paste(factor_cols, collapse = '+'))), List.test))
  }

  if (ncol(List.train) == K + 1) {
    xtrain = data.frame(x1 = List.train[, K + 1])
    xtest = data.frame(x1 = List.test[, K + 1])
  } else {
    xtrain = List.train[, -c(1:K)]
    xtest = List.test[, -c(1:K)]
  }

  xtrain = janitor::clean_names(xtrain)
  xtest = janitor::clean_names(xtest)

  if(num_cores == 1){
    fiti = tryCatch(SuperLearner(Y = as.numeric(List.train[,i]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitj = tryCatch(SuperLearner(Y = as.numeric(List.train[,j]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitij = tryCatch(SuperLearner(Y = as.numeric(pmin(List.train[, i], List.train[, j])),
                                      X = xtrain, family = binomial(), method = "method.AUC",
                                      SL.library = slib2, verbose = FALSE), silent = TRUE)

  }else{
    options(mc.cores = num_cores)
    getOption("mc.cores")
    cl <- parallel::makeCluster(num_cores, type = "PSOCK")
    parallel::clusterSetRNGStream(cl, iseed = 2343)
    foo <- parallel::clusterEvalQ(cl, library(SuperLearner))
    parallel::clusterExport(cl, foo)

    fiti = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(List.train[,i]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitj = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(List.train[,j]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitij = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(pmin(List.train[, i], List.train[, j])),
                                      X = xtrain, family = binomial(), method = "method.AUC",
                                      SL.library = slib2, verbose = FALSE), silent = TRUE)
    parallel::stopCluster(cl)
  }
  if ("try-error" %in% c(class(fiti), class(fitj), class(fitij))) {
    Warning("One or more fits with SuperLearner regression failed.")
    return(NULL)
  }
  else {
    q12 = pmax(pmin(predict(fitij, newdata = xtest, onlySL = TRUE)$pred,
                    1), eps)
    q1 = pmin(pmax(predict(fiti, newdata = xtest, onlySL = TRUE)$pred,
                   q12), 1)
    q2 = pmin(pmax(predict(fitj, newdata = xtest, onlySL = TRUE)$pred,
                   q12), 1)
    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
#' Estimate marginal and joint distribution of lists i and j using multinomial logistic model.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_mlogit(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_mlogit <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, ...){

  stopifnot(ncol(List.train) > K)

  require("nnet", quietly = TRUE)

  l = ncol(List.train) - K

  stopifnot(l > 0)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  eval(parse(text = paste0("List.train$caphist = paste0(List.train$L", i, ", List.train$L", j, ")")))
  eval(parse(text = paste0("List.test$caphist = paste0(List.test$L", i, ", List.test$L", j, ")")))

  mfit = try(multinom(formula = formula(paste("caphist ~ ", paste("x", 1:l, sep = '', collapse = ' + '))), data = List.train), silent = TRUE)
  if(!("try-error" %in% class(mfit))){

    pred = predict(mfit, newdata = List.test, "probs")
    if(is.null(setdiff(c("10", "11", "01"), colnames(pred)))){
      Error("Training data is missing one or more of the capture history combinations.")
    }
    q12 = pmax(pred[,"11"], eps)
    q1 = pmin(pred[,"10"] + q12, 1)
    q2 = pmin(pred[,"01"] + q12, 1)
  }else{
    Warning("One or more fits with SuperLearner regression failed.")
  }
  return(list(q1 = q1, q2 = q2, q12 = q12))
}
#' Estimate marginal and joint distribution of lists i and j using ensemble of ranger and logit.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param i The first list that is conditionally independent.
#' @param j The second list that is conditionally independent.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @return A list of the marginal and joint distribution probabilities q_1, q_2 and q_12.
#' @examples
#' qhat = qhat_ranger(List.train = List.train, List.test = List.test, K = 3, i = 1, j = 2, eps = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#'
#' @export
qhat_rangerlogit <- function(List.train, List.test, K = 2, i = 1, j = 2, eps = 0.005, ...){
  require("ranger", quietly = TRUE)
  l = ncol(List.train) - K

  stopifnot(l>0)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  fiti = ranger(formula(paste("factor(L", i, ") ~.", sep = '')), data = List.train[,-c(1:K)[-i]], probability = TRUE, classification = TRUE)
  fitj = ranger(formula(paste("factor(L", j, ") ~.", sep = '')), data = List.train[,-c(1:K)[-j]], probability = TRUE, classification = TRUE)
  fitij = ranger(formula(paste("factor(L", i, "*L", j, ") ~.", sep = '')), data = List.train[,c(i, j, K + 1:l)], probability = TRUE, classification = TRUE)

  if("try-error" %in% c(class(fiti), class(fitj), class(fitij))){
    Warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12.r = pmax(predict(fitij, data = List.test, type = "response")$predictions[,'1'], eps)
    q1.r = pmax(predict(fiti, data = List.test, type = "response")$predictions[,'1'], eps)
    q2.r = pmax(predict(fitj, data = List.test, type = "response")$predictions[,'1'], eps)

    ql = qhat_logit(List.train = List.train, List.test = List.test, K = K, i = i, j = j, eps = eps)
    q12.l = ql$q12
    q1.l = ql$q1
    q2.l = ql$q2

    require("nnls")
    A = cbind(1, q12.r, q12.l)
    coef = coef(nnls(b = List.test[,paste0('L',i)]*List.test[,paste0('L',j)], A = A))
    coef = coef/sum(coef)
    q12 = A %*% coef
    A = cbind(1, q1.r, q1.l)
    coef = coef(nnls(b = List.test[,paste0('L',i)], A = A))
    coef = coef/sum(coef)
    q1 = pmax(A %*% coef, q12)
    A = cbind(1, q2.r, q2.l)
    coef = coef(nnls(b = List.test[,paste0('L',j)], A = A))
    coef = coef/sum(coef)
    q2 = pmin(pmax(A %*% coef, q12), 1 - q1 + q12)

    #q12 = (lm(List.test[,paste0('L',i)]*List.test[,paste0('L',j)] ~ q12.r + q12.l))$fitted.values
    #q1 = (lm(List.test[,paste0('L',i)] ~ q1.r + q1.l))$fitted.values
    #q2 = (lm(List.test[,paste0('L',j)] ~ q2.r + q2.l))$fitted.values

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
