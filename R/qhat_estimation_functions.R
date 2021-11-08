#' Estimate marginal and joint distribution of lists j and k using logistic regression.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_logit(List.train = List.train, List.test = List.test, margin = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @importFrom stats binomial formula glm model.matrix na.omit predict rnorm sigma var
#' @export
qhat_logit <- function(List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005, ...){

  stopifnot(ncol(List.train) > K)
  if(missing(List.test)){
    List.test = List.train
  }
  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  fitj0 = try(glm(formula(paste("L", j, "*(1 - L", k, ") ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(j, k, (K + 1):ncol(List.train))]))
  fit0k = try(glm(formula(paste("L", k, "*(1 - L", j, ") ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(j, k, (K + 1):ncol(List.train))]))
  fitjk = try(glm(formula(paste("L", j, "*L", k, " ~.", sep = '')), family = binomial(link = "logit"), data = List.train[,c(j, k, (K + 1):ncol(List.train))]))

  if("try-error" %in% c(class(fitj0), class(fit0k), class(fitjk))){
    warning("One or more fits with logistic regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitjk, newdata = List.test, type = "response"), margin)
    q1 = pmax(q12 + predict(fitj0, newdata = List.test, type = "response"), margin)
    q2 = pmax(q12 + predict(fit0k, newdata = List.test, type = "response"), margin)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
#' Estimate marginal and joint distribution of lists j and k using generalized additive models.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_gam(List.train = List.train, List.test = List.test, margin = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @import gam
#' @references Trevor Hastie (2020). gam: Generalized Additive Models. _R package version 1.20_. https://CRAN.R-project.org/package=gam
#' @export
qhat_gam <- function(List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005, ...){

  stopifnot(ncol(List.train) > K)
  if(missing(List.test)){
    List.test = List.train
  }
  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  requireNamespace("gam", quietly = TRUE)
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

  fitj = try(gam::gam(formula(paste0("L", j, " ~", covformula)), data = List.train, family = binomial("logit")))
  fitk = try(gam::gam(formula(paste0("L", k, " ~", covformula)), data = List.train, family = binomial("logit")))
  fitjk = try(gam::gam(formula(paste0("L", j, "*L", k, " ~", covformula)), data = List.train, family = binomial("logit")))

  if("try-error" %in% c(class(fitj), class(fitk), class(fitjk))){
    warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitjk, newdata = List.test, type = "response"), margin)
    q1 = pmax(predict(fitj, newdata = List.test, type = "response"), margin)
    q2 = pmax(predict(fitk, newdata = List.test, type = "response"), margin)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
#' Estimate marginal and joint distribution of lists j and k using ranger.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_ranger(List.train = List.train, List.test = List.test, margin = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @import ranger
#' @references Marvin N. Wright, Andreas Ziegler (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. _Journal of Statistical Software_, *77*(1), 1-17. doi:10.18637/jss.v077.i01
#' @export
qhat_ranger <- function(List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005, ...){
  requireNamespace("ranger", quietly = TRUE)
  l = ncol(List.train) - K

  stopifnot(l>0)
  if(missing(List.test)){
    List.test = List.train
  }

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  fitj = ranger(formula(paste("factor(L", j, ") ~.", sep = '')), data = List.train[,-c(1:K)[-j]], probability = TRUE, classification = TRUE)
  fitk = ranger(formula(paste("factor(L", k, ") ~.", sep = '')), data = List.train[,-c(1:K)[-k]], probability = TRUE, classification = TRUE)
  fitjk = ranger(formula(paste("factor(L", j, "*L", k, ") ~.", sep = '')), data = List.train[,c(j, k, K + 1:l)], probability = TRUE, classification = TRUE)

  if("try-error" %in% c(class(fitj), class(fitk), class(fitjk))){
    warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12 = pmax(predict(fitjk, data = List.test, type = "response")$predictions[,'1'], margin)
    q1 = pmax(predict(fitj, data = List.test, type = "response")$predictions[,'1'], margin)
    q2 = pmax(predict(fitk, data = List.test, type = "response")$predictions[,'1'], margin)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
#' Estimate marginal and joint distribution of lists j and k using super learner.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param sl.lib The functions from the SuperLearner library to be used for model fitting. See [SuperLearner::listWrappers()].
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param num_cores The number of cores to be used for paralellization in Super Learner.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_sl(List.train = List.train, List.test = List.test, margin = 0.005, num_cores = 1)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' # One can specify the number of cores to be used for parallel computing
#' qhat = qhat_sl(List.train = List.train, List.test = List.test, margin = 0.005, num_cores = 2)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @import SuperLearner
#' @references Eric Polley, Erin LeDell, Chris Kennedy and Mark van der Laan (2021). SuperLearner: Super Learner Prediction. _R package version 2.0-28_. https://CRAN.R-project.org/package=SuperLearner
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, _Statistical Applications of Genetics and Molecular Biology_, 6, article 25.
#' @export
qhat_sl <- function (List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005,
                     sl.lib = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), num_cores = NA,...)
{

  stopifnot(ncol(List.train) > K)
  if(missing(List.test)){
    List.test = List.train
  }

  requireNamespace("SuperLearner", quietly = TRUE, warn.conflicts = FALSE)
  requireNamespace("parallel", quietly = TRUE, warn.conflicts = FALSE)
  requireNamespace("dplyr")
  requireNamespace("janitor", quietly = TRUE, warn.conflicts = FALSE)
  requireNamespace("tidyr", quietly = TRUE, warn.conflicts = FALSE)
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

  factor_cols <- subset(List.train, select = -c(1:K)) %>% select_if(Negate(is.numeric)) %>% names()

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
    fitj = tryCatch(SuperLearner(Y = as.numeric(List.train[,j]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitk = tryCatch(SuperLearner(Y = as.numeric(List.train[,k]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitjk = tryCatch(SuperLearner(Y = as.numeric(pmin(List.train[, j], List.train[, k])),
                                      X = xtrain, family = binomial(), method = "method.AUC",
                                      SL.library = slib2, verbose = FALSE), silent = TRUE)

  }else{
    options(mc.cores = num_cores)
    getOption("mc.cores")
    cl <- parallel::makeCluster(num_cores, type = "PSOCK")
    foo <-parallel::clusterEvalQ(cl, library(SuperLearner))
    parallel::clusterSetRNGStream(cl, iseed = 1)
    #
    parallel::clusterExport(cl, foo)

    fitj = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(List.train[,j]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitk = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(List.train[,k]), X = xtrain, family = binomial(),
                                     SL.library = slib2, method = "method.AUC",
                                     verbose = FALSE), silent = TRUE)
    fitjk = tryCatch(snowSuperLearner(cluster = cl, Y = as.numeric(pmin(List.train[, j], List.train[, k])),
                                      X = xtrain, family = binomial(), method = "method.AUC",
                                      SL.library = slib2, verbose = FALSE), silent = TRUE)
    parallel::stopCluster(cl)
  }
  if ("try-error" %in% c(class(fitj), class(fitk), class(fitjk))) {
    warning("One or more fits with SuperLearner regression failed.")
    return(NULL)
  }
  else {
    q12 = pmax(pmin(predict(fitjk, newdata = xtest, onlySL = TRUE)$pred,
                    1), margin)
    q1 = pmin(pmax(predict(fitj, newdata = xtest, onlySL = TRUE)$pred,
                   q12), 1)
    q2 = pmin(pmax(predict(fitk, newdata = xtest, onlySL = TRUE)$pred,
                   q12), 1)
    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
#' Estimate marginal and joint distribution of lists j and k using multinomial logistic model.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_mlogit(List.train = List.train, List.test = List.test, margin = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @importFrom nnet "multinom"
#' @references Croissant Y (2020). Estimation of Random Utility Models in R: The mlogit Package. _Journal of Statistical Software_, *95*(11), 1-41. doi: 10.18637/jss.v095.i11 (URL: https://doi.org/10.18637/jss.v095.i11).
#' @references Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. _Springer_, New York. ISBN 0-387-95457-0
#' @export
qhat_mlogit <- function(List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005, ...){

  stopifnot(ncol(List.train) > K)
  if(missing(List.test)){
    List.test = List.train
  }

  requireNamespace("nnet", quietly = TRUE)

  l = ncol(List.train) - K

  stopifnot(l > 0)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  eval(parse(text = paste0("List.train$caphist = paste0(List.train$L", j, ", List.train$L", k, ")")))
  eval(parse(text = paste0("List.test$caphist = paste0(List.test$L", j, ", List.test$L", k, ")")))

  mfit = try(multinom(formula = formula(paste("caphist ~ ", paste("x", 1:l, sep = '', collapse = ' + '))), data = List.train), silent = TRUE)
  if(!("try-error" %in% class(mfit))){

    pred = predict(mfit, newdata = List.test, "probs")
    if(is.null(setdiff(c("10", "11", "01"), colnames(pred)))){
      stop("Training data is missing one or more of the capture history combinations.")
    }
    q12 = pmax(pred[,"11"], margin)
    q1 = pmin(pred[,"10"] + q12, 1)
    q2 = pmin(pred[,"01"] + q12, 1)
  }else{
    warning("One or more fits with SuperLearner regression failed.")
  }
  return(list(q1 = q1, q2 = q2, q12 = q12))
}
#' Estimate marginal and joint distribution of lists j and k using ensemble of ranger and logit.
#'
#' @param List.train The training data matrix used to estimate the distibution functions.
#' @param List.test The data matrix on which the estimator function is applied.
#' @param K The number of lists in the data.
#' @param j The first list that is conditionally independent.
#' @param k The second list that is conditionally independent.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param ... Any extra arguments passed into the function.
#' @return A list of the marginal and joint distribution probabilities \code{q1}, \code{q2} and \code{q12}.
#' @examples
#' \dontrun{
#' qhat = qhat_ranger(List.train = List.train, List.test = List.test, margin = 0.005)
#' q1 = qhat$q1
#' q2 = qhat$q2
#' q12 = qhat$q12
#' }
#' @import ranger nnls
#' @references Marvin N. Wright, Andreas Ziegler (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. _Journal of Statistical Software_, *77*(1), 1-17. doi:10.18637/jss.v077.i01
#' @references Polley, Eric C. and van der Laan, Mark J., (May 2010) Super Learner In Prediction. _U.C. Berkeley Division of Biostatistics Working Paper Series_. Working Paper 266. https://biostats.bepress.com/ucbbiostat/paper266
#' @export
qhat_rangerlogit <- function(List.train, List.test, K = 2, j = 1, k = 2, margin = 0.005, ...){
  requireNamespace("ranger", quietly = TRUE)
  requireNamespace("nnls")
  l = ncol(List.train) - K
  if(missing(List.test)){
    List.test = List.train
  }

  stopifnot(l>0)

  colnames(List.train) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))
  colnames(List.test) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List.train) - K), sep = ''))

  fitj = ranger(formula(paste("factor(L", j, ") ~.", sep = '')), data = List.train[,-c(1:K)[-j]], probability = TRUE, classification = TRUE)
  fitk = ranger(formula(paste("factor(L", k, ") ~.", sep = '')), data = List.train[,-c(1:K)[-k]], probability = TRUE, classification = TRUE)
  fitjk = ranger(formula(paste("factor(L", j, "*L", k, ") ~.", sep = '')), data = List.train[,c(j, k, K + 1:l)], probability = TRUE, classification = TRUE)

  if("try-error" %in% c(class(fitj), class(fitk), class(fitjk))){
    warning("One or more fits with GAM regression failed.")
    return(NULL)
  }else{
    q12.r = pmax(predict(fitjk, data = List.test, type = "response")$predictions[,'1'], margin)
    q1.r = pmax(predict(fitj, data = List.test, type = "response")$predictions[,'1'], margin)
    q2.r = pmax(predict(fitk, data = List.test, type = "response")$predictions[,'1'], margin)

    ql = qhat_logit(List.train = List.train, List.test = List.test, K = K, j = j, k = k, margin = margin)
    q12.l = ql$q12
    q1.l = ql$q1
    q2.l = ql$q2

    requireNamespace("nnls")
    A = cbind(1, q12.r, q12.l)
    coef = coef(nnls(b = List.test[,paste0('L',j)]*List.test[,paste0('L',k)], A = A))
    coef = coef/sum(coef)
    q12 = A %*% coef
    A = cbind(1, q1.r, q1.l)
    coef = coef(nnls(b = List.test[,paste0('L',j)], A = A))
    coef = coef/sum(coef)
    q1 = pmax(A %*% coef, q12)
    A = cbind(1, q2.r, q2.l)
    coef = coef(nnls(b = List.test[,paste0('L',k)], A = A))
    coef = coef/sum(coef)
    q2 = pmin(pmax(A %*% coef, q12), 1 - q1 + q12)

    return(list(q1 = q1, q2 = q2, q12 = q12))
  }
}
