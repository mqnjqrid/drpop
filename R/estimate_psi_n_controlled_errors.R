#' Estimate the total population size and capture probabilities using perturbed true nuisance functions.
#'
#' @param List_matrix The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param n The true population size. Required to calculate the added error.
#' @param K The number of lists in the data. typically the first \code{K} rows of List_matrix.
#' @param sigma The standard deviation from zero of the added error.
#' @param alpha The rate of convergence.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param iter An integer denoting the maximum number of iterations allowed for targeted maximum likelihood method.
#' @param twolist The logical value of whether targeted maximum likelihood algorithm fits only two modes when K = 2.
#' @return A list of estimates containing the following components:
#' \item{psi}{  A matrix of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The rows represent the list pair which is assumed to be independent conditioned on the covariates.
#' The columns represent the model and method combinations (PI = plug-in, BC = bias-corrected, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
#' \item{sigma2}{  A matrix of the efficiency bound \code{sigma^2} in the same format as \code{psi}.}
#' \item{n}{  A matrix of the estimated population size n in the same format as \code{psi}.}
#' \item{varn}{  A matrix of the variance for population size estimate in the same format as \code{psi}.}
#' \item{N}{  The number of data points used in the estimation after removing rows with missing data.}
#'
#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2,1), nrow = nrow(data))
#'
#' psin_estimate = psinhat(List_matrix = data)
#' #this returns the basic plug-in estimate since covariates are absent.
#'
#' data = cbind(data, x)
#' psin_estimate = psinhat(List_matrix = data, funcname = c("logit", "sl"), nfolds = 2, twolist = FALSE, eps = 0.005)
#' #this returns the plug-in, the bias-corrected and the tmle estimate for the two models
#' @export
psinhat_simul = function(List_matrix, n, K, nfolds = 5, sigma, alpha, eps = 0.005, iter = 100, twolist = TRUE){

  if(missing(n)){
    n = nrow(List_matrix)
  }
  eta = n^alpha
  delta = 0
  #removing all rows with only 0's
  List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  List_matrix = as.data.frame(List_matrix)
  N = nrow(List_matrix)
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

  psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3)
  rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(k) {
    sapply((k + 1):K, function(s) {
      return(paste(k, ", ", s, sep = ''))
    })}))
  colnames(psiinv_summary) = c("PI", "BC", "TMLE")
  var_summary = psiinv_summary

  permutset = sample(1:N, N, replace = FALSE)

  for(i in 1:(K - 1)){
    if(!setequal(List_matrix[,i], c(0,1))){
      #     cat("List ", i, " is not in the required format or is degenerate.\n")
      next
    }
    for(j in (i + 1):K){
      if(!setequal(List_matrix[,j], c(0,1))){
        #       cat("List ", j, " is not in the required format or is degenerate.\n")
        next
      }
      psiinvmat = matrix(NA, nrow = nfolds, ncol = 3)
      colnames(psiinvmat) = c("PI", "BC", "TMLE")
      varmat = psiinvmat

      for(folds in 1:nfolds){

        if(nfolds == 1){
          List2 = as.data.frame(List_matrix)
        }else{
          sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
          sbset = sbset[sbset <= N]
          List2 = as.data.frame(List_matrix[permutset[sbset],])
        }

        xmat = as.matrix(List2[,-c(1:K)])
        p1 = unlist(apply(xmat, 1, pi1))
        p2 = unlist(apply(xmat, 1, pi2))

        q10_0 = p1*(1 - p2)/(1 - (1-p1)*(1-p2))
        q02_0 = p2*(1 - p1)/(1 - (1-p1)*(1-p2))
        q12_0 = p2*p1/(1 - (1-p1)*(1-p2))

        yi = List2[,paste("L", i, sep = '')]
        yj = List2[,paste("L", j, sep = '')]

        epsiln = matrix(rnorm(3*nrow(xmat), 1/eta, sigma/eta), ncol = 3)
        q12 = expit(logit(q12_0) + epsiln[,3])
        q1 = pmin(expit(logit(q10_0) + epsiln[,1]) + q12, 1)
        q2 = pmax(pmin(expit(logit(q02_0) + epsiln[,2]) + q12, 1 + q12 - q1), q12/q1)

        #colsubset = stringr::str_subset(colnames(psiinv_summary), func)

        gammainvhat = q1*q2/q12
        psiinvhat = mean(gammainvhat, na.rm = TRUE)

        phihat = gammainvhat*(yj/q2 + yi/q1 - yi*yj/q12) - psiinvhat

        Qnphihat = mean(phihat, na.rm = TRUE)

        psiinvhatq = max(psiinvhat + Qnphihat, 1)

        psiinvmat[folds, 1:2] = c(psiinvhat, psiinvhatq)

        sigmasq = var(phihat, na.rm = TRUE)
        varmat[folds, 1:2] = sigmasq/N

        datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12))
        datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
        colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")

        tmle = tmle(datmat = datmat, iter = iter, eps = eps, eps_stop = 0.00001, twolist = twolist)

        if(tmle$error){
          warning("TMLE did not run or converge.")
          psiinvmat[folds,colsubset][3] = NA
          varmat[folds,colsubset][3] = NA
        }else{
          datmat = tmle$datmat
          q12 = pmax(datmat$q12, eps)
          q1 = pmin(datmat$q12 + datmat$q10, 1)
          q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

          gammainvhat = q1*q2/q12
          psiinvhat = mean(gammainvhat, na.rm = TRUE)

          phihat = gammainvhat*(yi/q1 + yj/q2 - yi*yj/q12) - psiinvhat

          Qnphihat = mean(phihat, na.rm = TRUE)

          psiinvmat[folds, 3] = psiinvhat
          sigmasq = var(phihat, na.rm = TRUE)
          varmat[folds, 3] = sigmasq/N
        }
      }
      psiinv_summary[paste(i, ", ", j, sep = ''),] = colMeans(psiinvmat, na.rm = TRUE)
      var_summary[paste(i, ", ", j, sep = ''),] = colMeans(varmat, na.rm = TRUE)
    }
  }

  return(list(psi = 1/psiinv_summary, sigma2 = N*var_summary,
              n = N*psiinv_summary,
              varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1),
              N = N
  ))
}
