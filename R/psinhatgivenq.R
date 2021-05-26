#' Estimate total population size and capture probability using user provided set of models.
#'
#' @param List_matrix The data frame in capture-recapture format with two lists for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param qhateval A list object of class \code{qhateval} as returned by the \code{qhateval} function.
#' @param q1mat A dataframe with capture probabilities for the first list.
#' @param q2mat A dataframe with capture probabilities for the second list.
#' @param q12mat A dataframe with capture probabilities for both the lists simultaneously.
#' @param idfold The fold assignment of each row during estimation.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
#' @return A list of estimates containing the following components:
#' \item{psi}{  A dataframe of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The rows represent the list pair which is assumed to be independent conditioned on the covariates.
#' The columns represent the model and method combinations (PI = plug-in, DR = doubly-robust, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
#' \item{sigma2}{  A dataframe of the efficiency bound \code{sigma^2} in the same format as \code{psi}.}
#' \item{n}{  A dataframe of the estimated population size n in the same format as \code{psi}.}
#' \item{varn}{  A dataframe of the variance for population size estimate in the same format as \code{psi}.}
#' \item{N}{  The number of data points used in the estimation after removing rows with missing data.}
#' \item{ifvals}{  The estimated influence function values for the observed data. Each column corresponds to an element in funcname.}
#' \item{nuis}{  The estimated nuisance functions (q12, q1, q2) for each element in funcname.}
#' \item{nuistmle}{  The estimated nuisance functions (q12, q1, q2) from tmle for each element in funcname.}
#' \item{cin.l}{  The estimated lower bound of a 95% confidence interval of \code{n}.}
#' \item{cin.u}{  The estimated upper bound of a 95% confidence interval of \code{n}.}
#'
#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @examples
#' data = simuldata(1000, 1)$List_matrix
#' qhat_estimate = qhateval(List_matrix = data, funcname = c("logit", "gam"), nfolds = 2, eps = 0.005)
#' psin_estimate = psinhatgivenq(List_matrix = data, qhateval = qhat_estimate)
#' @export
psinhatgivenq <- function(List_matrix, i = 1, j = 2, eps = 0.005, qhateval, q1mat, q2mat, q12mat, idfold, TMLE = TRUE, ...){

  K = 2
  n = nrow(List_matrix)

  if(missing(i)){
    i = 1
  }
  if(missing(j)){
    j = 2
  }

  stopifnot(!(missing(qhateval) & missing(q1mat) & missing(q2mat) & missing(q12mat)))

  if(!missing(qhateval)){
    q1mat = qhateval$q1mat
    q2mat = qhateval$q2mat
    q12mat = qhateval$q12mat
    idfold = qhateval$idfold
  }

  stopifnot(!is.null(q1mat) & !is.null(q2mat) & !is.null(q12mat))

  funcname = colnames(q12mat)

  if(missing(idfold) | is.null(idfold)){
    idfold = rep(1, n)
  }
  nfolds = max(idfold)
  stopifnot(!is.null(dim(List_matrix)))

  if(!informat(List_matrix = List_matrix, K = K)){
    List_matrix <- reformat(List_matrix = List_matrix, capturelists = 1:K)
  }

  List_matrix = as.data.frame(List_matrix)
  #N = number of observed or captured units
  N = nrow(List_matrix)

  stopifnot(N > 1)

  #renaming the columns of List_matrix for ease of use
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

  psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3*length(funcname))
  rownames(psiinv_summary) = paste0(i, ",", j)
  colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  var_summary = psiinv_summary

  ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname))
  colnames(ifvals) = funcname
  rownames(ifvals) = rep(rownames(psiinv_summary), each = N)

  nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname))
  colnames(nuis) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
  rownames(nuis) = rownames(ifvals)
  nuistmle = nuis

  psiinvmat = matrix(NA, nrow = nfolds, ncol = 3*length(funcname))
  colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  varmat = psiinvmat

  for (folds in 1:nfolds){#print(folds)

    List2 = List_matrix[idfold == folds,]

    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]

    for (func in funcname){

      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)

      q12 = q12mat[idfold == folds, func]
      q1 = pmin(pmax(q12, q1mat[idfold == folds, func]), 1)
      q2 = pmax(q12/q1, pmin(q2mat[idfold == folds, func], 1 + q12 - q1, 1))

      nuis[idfold == folds, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

      gammainvhat = q1*q2/q12
      psiinvhat = mean(gammainvhat, na.rm = TRUE)

      phihat = gammainvhat*(yj/q2 + yi/q1 - yi*yj/q12) - psiinvhat
      ifvals[idfold == folds, func] = phihat

      Qnphihat = mean(phihat, na.rm = TRUE)

      psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

      psiinvmat[folds, paste(func, c("PI", "DR"), sep = '.')] = c(psiinvhat, psiinvhat.dr)

      sigmasq = var(phihat, na.rm = TRUE)
      varmat[folds, paste(func, c("PI", "DR"), sep = '.')] = sigmasq/N

      datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12))
      datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
      colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")

      if(TMLE){
        tmle = tmle(datmat = datmat, eps = eps, K = 2, ...)
      }else{
        tmle = list(error = TRUE)
      }

      if(tmle$error){
        warning("TMLE did not run or converge.")
        psiinvmat[folds, paste(func, "TMLE", sep = '.')] = NA
        varmat[folds, paste(func, "TMLE", sep = '.')] = NA
      }else{
        datmat = tmle$datmat
        q12 = pmax(datmat$q12, eps)
        q1 = pmin(datmat$q12 + datmat$q10, 1)
        q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

        nuistmle[idfold == folds, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

        gammainvhat = q1*q2/q12
        psiinvhat.tmle = mean(gammainvhat, na.rm = TRUE)

        phihat = gammainvhat*(yi/q1 + yj/q2 - yi*yj/q12) - psiinvhat.tmle

        Qnphihat = mean(phihat, na.rm = TRUE)

        psiinvmat[folds, paste(func, "TMLE", sep = '.')] = psiinvhat.tmle
        sigmasq = var(phihat, na.rm = TRUE)
        varmat[folds, paste(func, "TMLE", sep = '.')] = sigmasq/N
      }
    }
  }

  psiinv_summary[paste0(i, ",", j),] = colMeans(psiinvmat, na.rm = TRUE)
  var_summary[paste0(i, ",", j),] = colMeans(varmat, na.rm = TRUE)

  result <- list(psi = 1/psiinv_summary, sigma2 = N*var_summary, n = round(N*psiinv_summary),
              varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1), N = N,
              ifvals = ifvals, nuis = nuis, nuistmle = nuistmle,
              cin.l = round(pmax(N*psiinv_summary - 1.96*sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)), N)),
              cin.u = round(N*psiinv_summary + 1.96 *sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1))))
  class(result) <- "psinhat"
  return(result)
}
