#' Estimate total population size and capture probability using user provided set of models or user provided nuisance estimates.
#'
#' @param List_matrix The data frame in capture-recapture format with \code{K} lists for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the \code{K} lists. The remaining columns are covariates in numeric format.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param getnuis A list object with the nuisance function estimates and the fold assignment of the rows for cross-fitting.
#' @param q1mat A dataframe with capture probabilities for the first list.
#' @param q2mat A dataframe with capture probabilities for the second list.
#' @param q12mat A dataframe with capture probabilities for both the lists simultaneously.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param idfold The fold assignment of each row during estimation.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
#' @param PLUGIN The logical value to indicate whether the plug-in estimates is returned.
#' @return A list of estimates containing the following components for each list-pair, model and method (PI = plug-in, DR = doubly-robust, TMLE = targeted maximum likelihood estimate):
#' \item{result}{  A dataframe of the below estimated quantities.
#' \itemize{
#' \item{psi}{  The estimated capture probability.}
#' \item{sigma}{  The efficiency bound.}
#' \item{n}{  The estimated population size n.}
#' \item{sigman}{  The estimated standard deviation of the population size.}
#' \item{cin.l}{  The estimated lower bound of a 95% confidence interval of \code{n}.}
#' \item{cin.u}{  The estimated upper bound of a 95% confidence interval of \code{n}.}}}
#' \item{N}{  The number of data points used in the estimation after removing rows with missing data.}
#' \item{ifvals}{  The estimated influence function values for the observed data.}
#' \item{nuis}{  The estimated nuisance functions (q12, q1, q2) for each element in funcname.}
#' \item{nuistmle}{  The estimated nuisance functions (q12, q1, q2) from tmle for each element in funcname.}
#' \item{idfold}{  The division of the rows into sets (folds) for cross-fitting.}
#'
#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @examples
#' data = simuldata(1000, 1)$List_matrix
#' qhat_estimate = getnuis(List_matrix = data, funcname = c("logit", "gam"), nfolds = 2, eps = 0.005)
#' psin_estimate = popsize(List_matrix = data, getnuis = qhat_estimate)
#' @export
popsize <- function(List_matrix, j = 1, k = 2, eps = 0.005, getnuis, q1mat, q2mat, q12mat, idfold, TMLE = TRUE, PLUGIN = TRUE, ...){

  if(missing(getnuis) & missing(q1mat) & missing(q2mat) & missing(q12mat)){
    return(popsize_base(List_matrix, eps = eps, TMLE = TMLE, PLUGIN = PLUGIN, ...))
  }

  K = 2
  n = nrow(List_matrix)

  if(missing(j)){
    j = 1
  }
  if(missing(k)){
    k = 2
  }

  if(!missing(getnuis)){
    q1mat = getnuis$q1mat
    q2mat = getnuis$q2mat
    q12mat = getnuis$q12mat
    idfold = getnuis$idfold
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
  rownames(psiinv_summary) = paste0(j, ",", k)
  colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  var_summary = psiinv_summary

  ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname) + 1)
  colnames(ifvals) = c("listpair", funcname)
  ifvals[,"listpair"] = rep(rownames(psiinv_summary), each = N)

  nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname) + 1)
  colnames(nuis) = c("listpair", paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.'))
  nuis = as.data.frame(nuis)
  sapply(nuis, "class")
  nuis[,"listpair"] = ifvals[,"listpair"]
  if(TMLE){
    nuistmle = nuis
  }

  psiinvmat = matrix(NA, nrow = nfolds, ncol = 3*length(funcname))
  colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  varmat = psiinvmat

  for (folds in 1:nfolds){#print(folds)

    List2 = List_matrix[idfold == folds,]

    yj = List2[,paste("L", j, sep = '')]
    yk = List2[,paste("L", k, sep = '')]

    for (func in funcname){

      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)

      q12 = q12mat[idfold == folds, func]
      q1 = pmin(pmax(q12, q1mat[idfold == folds, func]), 1)
      q2 = pmax(q12/q1, pmin(q2mat[idfold == folds, func], 1 + q12 - q1, 1))

      nuis[idfold == folds, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

      gammainvhat = q1*q2/q12
      psiinvhat = mean(gammainvhat, na.rm = TRUE)

      phihat = gammainvhat*(yk/q2 + yj/q1 - yj*yk/q12) - psiinvhat
      ifvals[idfold == folds, func] = phihat

      Qnphihat = mean(phihat, na.rm = TRUE)

      psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

      psiinvmat[folds, paste(func, c("PI", "DR"), sep = '.')] = c(psiinvhat, psiinvhat.dr)

      sigmasq = var(phihat, na.rm = TRUE)
      varmat[folds, paste(func, c("PI", "DR"), sep = '.')] = sigmasq/N

      datmat = as.data.frame(cbind(yj, yk, yj*yk, q1 - q12, q2 - q12, q12))
      datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
      colnames(datmat) = c("yj", "yk", "yjk", "q10", "q02", "q12")

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

        phihat = gammainvhat*(yj/q1 + yk/q2 - yj*yk/q12) - psiinvhat.tmle

        Qnphihat = mean(phihat, na.rm = TRUE)

        psiinvmat[folds, paste(func, "TMLE", sep = '.')] = psiinvhat.tmle
        sigmasq = var(phihat, na.rm = TRUE)
        varmat[folds, paste(func, "TMLE", sep = '.')] = sigmasq/N
      }
    }
  }

  psiinv_summary[paste0(j, ",", k),] = colMeans(psiinvmat, na.rm = TRUE)
  var_summary[paste0(j, ",", k),] = colMeans(varmat, na.rm = TRUE)

  result <- list(psi = 1/psiinv_summary, sigma = sqrt(N*var_summary), n = round(N*psiinv_summary),
                 sigman = sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)),
                 cin.l = round(pmax(N*psiinv_summary - 1.96*sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)), N)),
                 cin.u = round(N*psiinv_summary + 1.96 *sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1))))
  result <- Reduce(function(...) merge(..., by = c("listpair", "Var2")),
                   lapply(1:length(result), function(i)
                     reshape2::melt(result[[i]], value.name = names(result)[i], varnames = c("listpair", "Var2"))))
  result <- tidyr::separate(data = result, col = "Var2", into = c("model", "method"), sep = '\\.')
  if(!TMLE){
    result = result[result$method != "TMLE",]
  }
  if(!PLUGIN){
    result = result[result$method != "PI",]
  }else{
    warning("Plug-in variance is not well-defined. Returning variance evaluated using DR estimator formula")
  }
  ifvals = as.data.frame(ifvals)
  ifvals$listpair = paste0(j, ',', k)
  nuis = as.data.frame(nuis)
  nuis$listpair = paste0(j, ',', k)

  object = list(result = result, N = N, ifvals = as.data.frame(ifvals), nuis = as.data.frame(nuis), idfold = idfold)
  if(TMLE){
    nuistmle = as.data.frame(nuistmle)
    nuistmle$listpair = paste0(j, ',', k)
    object$nuistmle = as.data.frame(nuistmle)
  }
  class(object) = "popsize"
  return(invisible(object))
}
