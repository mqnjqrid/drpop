#' Estimate total population size and capture probability using user provided set of models or user provided nuisance estimates.
#'
#' @param data The data frame in capture-recapture format with \code{K} lists for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the \code{K} lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists that are present in the data.
#' @param j The first list to be used for estimation.
#' @param k The secod list to be used in the estimation.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param getnuis A list object with the nuisance function estimates and the fold assignment of the rows for cross-fitting or a data.frame with the nuisance estimates.
#' @param q1mat A dataframe with capture probabilities for the first list.
#' @param q2mat A dataframe with capture probabilities for the second list.
#' @param q12mat A dataframe with capture probabilities for both the lists simultaneously.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param sl.lib algorithm library for SuperLearner. Default library includes "gam", "glm", "glmnet", "glm.interaction", "ranger".
# (See \code{\link[pkg:SuperLearner]{listWrappers}})
#' @param Nmin The cutoff for minimum sample size to perform doubly robust estimation. Otherwise, Petersen estimator is returned.
#' @param idfold The fold assignment of each row during estimation.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
#' @param PLUGIN The logical value to indicate whether the plug-in estimates is returned.
#' @param ... Any extra arguments passed into the function.
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
#' @references Bickel, P. J., Klaassen, C. A., Bickel, P. J., Ritov, Y., Klaassen, J., Wellner, J. A., and Ritov, Y. (1993). Efficient and adaptive estimation for semiparametric models, volume 4. _Johns Hopkins University Press Baltimore_
#' @references van der Vaart, A. (2002a). Part iii: Semiparameric statistics. Lectures on Probability Theory and Statistics, pages 331-457
#' @references van der Laan, M. J. and Robins, J. M. (2003). Unified methods for censored longitudinal data and causality. _Springer Science & Business Media_
#' @references Tsiatis, A. (2006). Semiparametric theory and missing data _springer. New York_
#' @references Kennedy, E. H. (2016). Semiparametric theory and empirical processes in causal inference. _Statistical causal inferences and their applications in public health research_, pages 141-167. _Springer_
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust capture-recapture methods for estimating population size. _arXiv preprint_ *arXiv:2104.14091*.
#' @examples
#' \donttest{
#' data = simuldata(1000, l = 3)$data
#' qhat = popsize(data = data, funcname = c("logit", "gam"), nfolds = 2, margin = 0.005)
#' psin_estimate = popsize(data = data, getnuis = qhat$nuis, idfold = qhat$idfold)
#'
#' data = simuldata(n = 6000, l = 3)$data
#' psin_estimate = popsize(data = data[,1:2])
#' #this returns the basic plug-in estimate since covariates are absent.
#'
#' psin_estimate = popsize(data = data, funcname = c("gam", "rangerlogit"))
#' }
#' @importFrom dplyr "%>%" "mutate" "select_if"
#' @import utils
#' @import stats
#' @importFrom tidyr "separate"
#' @export
popsize <- function(data, K = 2, j, k, margin = 0.005, filterrows = FALSE, nfolds = 5, funcname = c("rangerlogit"),
                     sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), getnuis, q1mat, q2mat, q12mat, idfold, TMLE = TRUE, PLUGIN = TRUE, Nmin = 100,...){

  if(!missing(j) & !missing(k)){
    if(j == k) {
      k = j %% K + 1
      warning(paste0("Selected lists are identical. Using k = ", k, "."))
    }
    if(j > k){
      j = j + k
      k = j - k
      j = j - k
      warning("Switching j and k to ensure j < k.")
    }
  }

  if(missing(getnuis) & missing(q1mat) & missing(q2mat) & missing(q12mat)){
    if(!missing(j) & missing(k)){
      if(j == K)
        return(popsize_base(data, K = K, k0 = j, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                          sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...))
      else
        return(popsize_base(data, K = K, j0 = j, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                            sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...))
    }else if(missing(j) & !missing(k)){
      if(k < K)
        return(popsize_base(data, K = K, j0 = k, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                            sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...))
      else
        return(popsize_base(data, K = K, k0 = k, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                            sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...))
    }else
      return(popsize_base(data, K = K, j0 = j, k0 = k, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                        sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...))
  }

  K = 2
  n = nrow(data)

  if(missing(j)){
    j = 1
  }
  if(missing(k)){
    k = 2
  }

  if(!missing(getnuis)){
    if(class(getnuis) == "data.frame"){
      q1mat = subset(getnuis, select = grep(colnames(getnuis), pattern = 'q1$'))
      colnames(q1mat) = stringr::str_remove_all(colnames(q1mat), '\\.q1')
      q2mat = subset(getnuis, select = grep(colnames(getnuis), pattern = 'q2$'))
      colnames(q2mat) = stringr::str_remove_all(colnames(q1mat), '\\.q2')
      q12mat = subset(getnuis, select = grep(colnames(getnuis), pattern = 'q12$'))
      colnames(q12mat) = stringr::str_remove_all(colnames(q1mat), '\\.q12')
    }else{
      q1mat = getnuis$q1mat
      q2mat = getnuis$q2mat
      q12mat = getnuis$q12mat
      idfold = getnuis$idfold
    }
  }

  stopifnot(!is.null(q1mat) & !is.null(q2mat) & !is.null(q12mat))

  funcname = colnames(q12mat)

  if(missing(idfold) | is.null(idfold)){
    idfold = rep(1, n)
  }
  nfolds = max(idfold)
  stopifnot(!is.null(dim(data)))

  if(!informat(data = data, K = K)){
    data <- reformat(data = data, capturelists = 1:K)
  }

  data = as.data.frame(data)
  #N = number of observed or captured units
  N = nrow(data)

  stopifnot(N > 1)

  #renaming the columns of data for ease of use
  colnames(data) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(data) - K), sep = ''))

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

  for (folds in 1:nfolds){

    List2 = data[idfold == folds,]

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
      datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, margin), 1 - margin))}))
      colnames(datmat) = c("yj", "yk", "yjk", "q10", "q02", "q12")

      if(TMLE){
        tmle = tmle(datmat = datmat, margin = margin, K = 2, ...)
      }else{
        tmle = list(error = TRUE)
      }

      if(tmle$error){
        warning("TMLE did not run or converge.")
        psiinvmat[folds, paste(func, "TMLE", sep = '.')] = NA
        varmat[folds, paste(func, "TMLE", sep = '.')] = NA
      }else{
        datmat = tmle$datmat
        q12 = pmax(datmat$q12, margin)
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

popsize_base <- function(data, K = 2, j0, k0, filterrows = FALSE, funcname = c("rangerlogit"), nfolds = 5, margin = 0.005,
                         sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), Nmin = 500, TMLE = TRUE, PLUGIN = TRUE,...){

  requireNamespace("dplyr", quietly = TRUE, warn.conflicts = FALSE)
  requireNamespace("tidyr")
  l = ncol(data) - K
  n = nrow(data)

  stopifnot(!is.null(dim(data)))

  if(!informat(data = data, K = K)){
    data <- reformat(data = data, capturelists = 1:K)
  }

  data = na.omit(data)

  if(filterrows){
    #removing all rows with only 0's
    data = data[which(rowSums(data[,1:K]) > 0),]
  }

  data = as.data.frame(data)
  #N = number of observed or captured units
  N = nrow(data)

  stopifnot(N > 1)

  if (l >= 0 & N < Nmin){
    l = 0
    warning(cat("Insufficient number of observations for doubly-robust estimation."))
  }

  conforminglists = apply(data[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    stop("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }

  if(sum(conforminglists) < K){
    message(cat("Lists ", which(conforminglists == FALSE), " are not in the required format."))
  }

  if(!missing(j0))
    list1_vec = j0
  else
    list1_vec = c(1:(K-1))

  if(!missing(k0))
    list2_vec = k0
  else
    list2_vec = c(1:K)

  if(l == 0){
    #renaming the columns of data for ease of use
    colnames(data) = c(paste("L", 1:K, sep = ''))

    listpair = unlist(sapply(list1_vec, function(j1) {
      sapply(setdiff(list2_vec, list1_vec[list1_vec <= j1]), function(k1) {
        return(paste(min(j1, k1), ",", max(j1, k1), sep = ''))
      })}))
    psiinv = data.frame(listpair = listpair)
    psiinv$psiin = NA
    psiinv$sigma = NA

    for(j in list1_vec){
      j0 = j
      if(!setequal(data[,j], c(0,1))){
        next
      }
      for(k in setdiff(list2_vec, list1_vec[list1_vec <= j0])){
        if(!setequal(data[,k], c(0,1))){
          next
        }
        if(j0 > k) {
          j = k
          k = j0
        }else
          j = j0
        q1 = mean(data[,j])
        q2 = mean(data[,k])
        q12 = mean(data[,j]*data[,k])
        psiinv[psiinv$listpair == paste0(j, ",", k),]$psiin = pmax(q1*q2/q12, 1)
        psiinv[psiinv$listpair == paste0(j, ",", k),]$sigma = sqrt(q1*q2*pmax(q1*q2 - q12, 0)*(1 - q12)/q12^3/N)
      }
    }

    result <- psiinv %>% mutate(psi = 1/psiin, sigma = sqrt(N)*sigma, n = round(N*psiin),
                                sigman = sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)),
                                cin.l = round(pmax(N*psiin - 1.96*sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)), N)),
                                cin.u = round(N*psiin + 1.96 *sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)))) %>% as.data.frame()

    result = subset(result, select = -c(psiin))

    object = list(result = result, N = N)
    class(object) = "popsize"
    return(object)
  }else{

    #renaming the columns of data for ease of use
    colnames(data) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(data) - K), sep = ''))

    if(nfolds > 1 & nfolds > N/50) {
      nfolds = pmax(floor(N/50), 1)
      cat("nfolds is reduced to ", nfolds, " to have sufficient test data.\n")
    }

    listpair = unlist(sapply(list1_vec, function(j1) {
      sapply(setdiff(list2_vec, list1_vec[list1_vec <= j1]), function(k1) {
        return(paste(min(j1, k1), ",", max(j1, k1), sep = ''))
      })}))
    psiinv_summary = matrix(0, nrow = length(listpair), ncol = 3*length(funcname))
    rownames(psiinv_summary) = listpair
    colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
    var_summary = psiinv_summary

    ifvals = matrix(NA, nrow = N*length(listpair), ncol = length(funcname) + 1)
    colnames(ifvals) = c("listpair", funcname)
    ifvals[,"listpair"] = rep(rownames(psiinv_summary), each = N)

    nuis = matrix(NA, nrow = N*length(listpair), ncol = 3*length(funcname) + 1)
    colnames(nuis) = c("listpair", paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.'))
    nuis = as.data.frame(nuis)
    sapply(nuis, "class")
    nuis[,"listpair"] = ifvals[,"listpair"]

    if(TMLE){
      nuistmle = nuis
    }

    permutset = sample(1:N, N, replace = FALSE)

    for(j in list1_vec){
      j0 = j
      if(!setequal(data[,j], c(0,1))){
        #     cat("List ", j, " is not in the required format or is degenerate.\n")
        next
      }
      for(k in setdiff(list2_vec, list1_vec[list1_vec <= j0])){
        if(!setequal(data[,k], c(0,1))){
          #       cat("List ", k, " is not in the required format or is degenerate.\n")
          next
        }
        if(j0 > k){
          j = k
          k = j0
        }else
          j = j0
        psiinvmat = matrix(numeric(0), nrow = nfolds, ncol = 3*length(funcname))
        colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
        varmat = psiinvmat

        ifvalsfold = matrix(numeric(0), nrow = N, ncol = length(funcname))
        colnames(ifvalsfold) = funcname

        nuisfold = matrix(numeric(0), nrow = N, ncol = 3*length(funcname))
        colnames(nuisfold) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
        nuistmlefold = nuisfold

        idfold = rep(1, N)

        for(folds in 1:nfolds){

          if(nfolds == 1){
            List1 = data
            List2 = List1
            sbset = 1:N
          }else{
            sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
            sbset = sbset[sbset <= N]
            List1 = data[permutset[-sbset],]
            List2 = data[permutset[sbset],]
            idfold[permutset[sbset]] = folds
          }

          yj = List2[,paste("L", j, sep = '')]
          yk = List2[,paste("L", k, sep = '')]

          overlapjk = mean(List1[,j]*List1[,k])
          if(overlapjk < margin) {
            warning(cat("Overlap between the lists ", j, " and ", k, " is less than ", margin, '.\n', sep = ''))
          }
          for (func in funcname){

            #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
            qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, j, k, margin = margin,...), silent = TRUE)

            if ("try-error" %in% class(qhat)) {
              next
            }

            q12 = qhat$q12
            q1 = pmin(pmax(q12, qhat$q1), 1)
            q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))

            nuisfold[permutset[sbset], paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

            gammainvhat = q1*q2/q12
            psiinvhat = mean(gammainvhat, na.rm = TRUE)

            phihat = gammainvhat*(yk/q2 + yj/q1 - yj*yk/q12) - psiinvhat
            ifvalsfold[permutset[sbset], func] = phihat

            Qnphihat = mean(phihat, na.rm = TRUE)

            psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

            psiinvmat[folds, paste(func, c("PI", "DR"), sep = '.')] = c(psiinvhat, psiinvhat.dr)

            sigmasq = var(phihat, na.rm = TRUE)
            varmat[folds, paste(func, c("PI", "DR"), sep = '.')] = sigmasq/N

            datmat = as.data.frame(cbind(yj, yk, yj*yk, q1 - q12, q2 - q12, q12))
            colnames(datmat) = c("yj", "yk", "yjk", "q10", "q02", "q12")

            if(TMLE) {
              tmle = tmle(datmat = datmat, margin = margin, K = K,...)
            }else{
              next
            }

            if(tmle$error){
              warning("TMLE did not run or converge.")
              psiinvmat[folds, paste(func, "TMLE", sep = '.')] = NA
              varmat[folds, paste(func, "TMLE", sep = '.')] = NA
            }else{
              datmat = tmle$datmat
              q12 = pmax(datmat$q12, margin)
              q1 = pmin(datmat$q12 + datmat$q10, 1)
              q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

              nuistmlefold[permutset[sbset], paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

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

        ifvals[ifvals[,"listpair"] == paste0(j, ",", k), colnames(ifvals) != "listpair"] = ifvalsfold
        nuis[nuis[, "listpair"] == paste0(j, ",", k), colnames(nuis) != "listpair"] = nuisfold
        if(TMLE){
          nuistmle[nuistmle[, "listpair"] == paste0(j, ",", k), colnames(nuistmle) != "listpair"] = nuistmlefold
        }
      }
    }

    result <- list(psi = 1/psiinv_summary, sigma = sqrt(N*var_summary), n = round(N*psiinv_summary),
                   sigman = sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)),
                   cin.l = round(pmax(N*psiinv_summary - 1.96*sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)), N)),
                   cin.u = round(N*psiinv_summary + 1.96 *sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1))))
    result <- Reduce(function(...) merge(..., by = c("listpair", "Var2")),
                     lapply(1:length(result), function(i)
                       reshape2::melt(result[[i]], value.name = names(result)[i], varnames = c("listpair", "Var2"))))
    result <- separate(data = result, col = "Var2", into = c("model", "method"), sep = '\\.')
    if(!TMLE){
      result = result[result$method != "TMLE",]
    }
    if(!PLUGIN){
      result = result[result$method != "PI",]
    }else{
      warning("Plug-in variance is not well-defined. Returning variance evaluated using DR estimator formula")
    }
    object = list(result = result, N = N, ifvals = as.data.frame(ifvals), nuis = as.data.frame(nuis), idfold = idfold)
    if(TMLE){
      object$nuistmle = as.data.frame(nuistmle)
    }
    class(object) = "popsize"
    return(object)
  }
}
#' @export
print.popsize <- function(x, ...){
  x$result$psi = round(x$result$psi, 3)
  x$result$sigma = round(x$result$sigma, 3)
  x$result$sigman = round(x$result$sigman, 3)
  print(x$result)
  invisible(x)
}
