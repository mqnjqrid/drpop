#' Estimate the total population size and capture probabilities using perturbed true nuisance functions.
#'
#' @param data The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param n The true population size. Required to calculate the added error.
#' @param K The number of lists in the data. typically the first \code{K} rows of data.
#' @param omega The standard deviation from zero of the added error.
#' @param alpha The rate of convergence. Takes values in (0, 1].
#' @param nfolds The number of folds to be used for cross fitting.
#' @param pi1 The function to calculate the conditional capture probabilities of list 1 using covariates.
#' @param pi2 The function to calculate the conditional capture probabilities of list 2 using covariates.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param iter An integer denoting the maximum number of iterations allowed for targeted maximum likelihood method.
#' @param twolist The logical value of whether targeted maximum likelihood algorithm fits only two modes when K = 2.
#' @return A list of estimates containing the following components:
#' \item{psi}{  A matrix of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The rows represent the list pair which is assumed to be independent conditioned on the covariates.
#' The columns represent the model and method combinations (PI = plug-in, DR = bias-corrected, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
#' \item{sigma2}{  A matrix of the efficiency bound \code{sigma^2} in the same format as \code{psi}.}
#' \item{n}{  A matrix of the estimated population size n in the same format as \code{psi}.}
#' \item{varn}{  A matrix of the variance for population size estimate in the same format as \code{psi}.}
#' \item{N}{  The number of data points used in the estimation after removing rows with missing data.}
#'
#' @examples
#' simulresult = simuldata(n = 2000, l = 2)
#' data = simulresult$data
#'
#' psin_estimate = popsize_simul(data = data,
#'       pi1 = simulresult$pi1, pi2 = simulresult$pi2,
#'       alpha = 0.25, omega = 1)
#'
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust capture-recapture methods for estimating population size. _arXiv preprint_ *arXiv:2104.14091*
#' @export
popsize_simul = function(data, n, K = 2, nfolds = 5, pi1, pi2, omega, alpha, margin = 0.005, iter = 100, twolist = TRUE){

  if(missing(n)){
    n = nrow(data)
  }
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  logit = function(x) {
    log(x/(1 - x))
  }
  eta = n^alpha
  delta = 0
  #removing all rows with only 0's
  data = data[which(rowSums(data[,1:K]) > 0),]
  data = as.data.frame(data)
  N = nrow(data)
  colnames(data) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(data) - K), sep = ''))

  psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3)
  rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(r) {
    sapply((r + 1):K, function(s) {
      return(paste(r, ", ", s, sep = ''))
    })}))
  colnames(psiinv_summary) = c("PI", "DR", "TMLE")
  var_summary = psiinv_summary

  permutset = sample(1:N, N, replace = FALSE)

  for(j in 1:(K - 1)){
    if(!setequal(data[,j], c(0,1))){
      #     cat("List ", j, " is not in the required format or is degenerate.\n")
      next
    }
    for(k in (j + 1):K){
      if(!setequal(data[,k], c(0,1))){
        #       cat("List ", k, " is not in the required format or is degenerate.\n")
        next
      }
      psiinvmat = matrix(NA, nrow = nfolds, ncol = 3)
      colnames(psiinvmat) = c("PI", "DR", "TMLE")
      varmat = psiinvmat

      for(folds in 1:nfolds){

        if(nfolds == 1){
          List2 = as.data.frame(data)
        }else{
          sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
          sbset = sbset[sbset <= N]
          List2 = as.data.frame(data[permutset[sbset],])
        }

        xmat = as.matrix(List2[,-c(1:K)])
        p1 = unlist(sapply(rowSums(xmat), pi1))
        p2 = unlist(sapply(rowSums(xmat), pi2))

        q10_0 = p1*(1 - p2)/(1 - (1-p1)*(1-p2))
        q02_0 = p2*(1 - p1)/(1 - (1-p1)*(1-p2))
        q12_0 = p2*p1/(1 - (1-p1)*(1-p2))

        yj = List2[,paste("L", j, sep = '')]
        yk = List2[,paste("L", k, sep = '')]

        marginiln = matrix(rnorm(3*nrow(xmat), 1/eta, omega/eta), ncol = 3)
        q12 = expit(logit(q12_0) + marginiln[,3])
        q1 = pmin(expit(logit(q10_0) + marginiln[,1]) + q12, 1)
        q2 = pmax(pmin(expit(logit(q02_0) + marginiln[,2]) + q12, 1 + q12 - q1), q12/q1)

        #colsubset = stringr::str_subset(colnames(psiinv_summary), func)

        gammainvhat = q1*q2/q12
        psiinvhat = mean(gammainvhat, na.rm = TRUE)

        phihat = gammainvhat*(yk/q2 + yj/q1 - yj*yk/q12) - psiinvhat

        Qnphihat = mean(phihat, na.rm = TRUE)

        psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

        psiinvmat[folds, 1:2] = c(psiinvhat, psiinvhat.dr)

        sigmasq = var(phihat, na.rm = TRUE)
        varmat[folds, 1:2] = sigmasq/N

        datmat = as.data.frame(cbind(yj, yk, yj*yk, q1 - q12, q2 - q12, q12))
        datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, margin), 1 - margin))}))
        colnames(datmat) = c("yj", "yk", "yjk", "q10", "q02", "q12")

        tmle = tmle(datmat = datmat, iter = iter, margin = margin, margin_stop = 0.00001, twolist = twolist)

        if(tmle$error){
          warning("TMLE did not run or converge.")
          psiinvmat[folds,3] = NA
          varmat[folds,3] = NA
        }else{
          datmat = tmle$datmat
          q12 = pmax(datmat$q12, margin)
          q1 = pmin(datmat$q12 + datmat$q10, 1)
          q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

          gammainvhat = q1*q2/q12
          psiinvhat.dr = mean(gammainvhat, na.rm = TRUE)

          phihat = gammainvhat*(yj/q1 + yk/q2 - yj*yk/q12) - psiinvhat.dr

          Qnphihat = mean(phihat, na.rm = TRUE)

          psiinvmat[folds, 3] = psiinvhat.dr
          sigmasq = var(phihat, na.rm = TRUE)
          varmat[folds, 3] = sigmasq/N
        }
      }
      psiinv_summary[paste(j, ", ", k, sep = ''),] = colMeans(psiinvmat, na.rm = TRUE)
      var_summary[paste(j, ", ", k, sep = ''),] = colMeans(varmat, na.rm = TRUE)
    }
  }

  return(list(psi = 1/psiinv_summary, sigma2 = N*var_summary,
              n = N*psiinv_summary,
              varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1),
              N = N
  ))
}
