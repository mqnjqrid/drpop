#' Estimate total population size and capture probability using user provided set of models.
#'
#' @param List_matrix The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists in the data. typically the first \code{K} rows of List_matrix.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param sl.lib algorithm library for SuperLearner. Default library includes "gam", "glm", "glmnet", "glm.interaction", "ranger".
#' @param Nmin The cutoff for minimum sample size to perform doubly robust estimation. Otherwise, Petersen estimator is returned.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
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
#setClass("psinhat", contains = "list")
# @exportClass psinhat
#setMethod("print", "psinhat", print.psinhat)
#' @export
psinhat <- function(List_matrix, K = 2, filterrows = FALSE, funcname = c("rangerlogit"), nfolds = 5, eps = 0.005,
                    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), Nmin = 500, TMLE = TRUE, ...){

  require("tidyverse", quietly = TRUE, warn.conflicts = FALSE)
  l = ncol(List_matrix) - K
  n = nrow(List_matrix)

  stopifnot(!is.null(dim(List_matrix)))

  if(!informat(List_matrix = List_matrix, K = K)){
    List_matrix <- reformat(List_matrix = List_matrix, capturelists = 1:K)
   }

  List_matrix = na.omit(List_matrix)

  if(filterrows){
    #removing all rows with only 0's
    List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  }

  List_matrix = as.data.frame(List_matrix)
  #N = number of observed or captured units
  N = nrow(List_matrix)

  stopifnot(N > 1)

  if (l >= 0 & N < Nmin){
    l = 0
    warning(cat("Insufficient number of observations for doubly-robust estimation."))
  }

  conforminglists = apply(List_matrix[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    stop("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }

  if(sum(conforminglists) < K){
    Message(cat("Lists ", which(conforminglists == FALSE), " are not in the required format."))
  }

  if(l == 0){
    #renaming the columns of List_matrix for ease of use
    colnames(List_matrix) = c(paste("L", 1:K, sep = ''))

    listpair = unlist(sapply(1:(K - 1), function(k) {
      sapply((k + 1):K, function(s) {
        return(paste(k, ",", s, sep = ''))
      })}))
    psiinv = data.frame(listpair = listpair)
    psiinv$psiin = NA
    psiinv$sigma = NA

    for(i in 1:(K - 1)){
      if(!setequal(List_matrix[,i], c(0,1))){
        next
      }
      for(j in (i + 1):K){
        if(!setequal(List_matrix[,j], c(0,1))){
          next
        }
        q1 = mean(List_matrix[,i])
        q2 = mean(List_matrix[,j])
        q12 = mean(List_matrix[,i]*List_matrix[,j])
        psiinv[psiinv$listpair == paste0(i, ",", j),]$psiin = q1*q2/q12
        psiinv[psiinv$listpair == paste0(i, ",", j),]$sigma = sqrt(q1*q2*(q1*q2 - q12)*(1 - q12)/q12^3/N)
      }
    }

    result <- psiinv %>% mutate(psi = 1/psiin, sigma = sqrt(N)*sigma, n = round(N*psiin),
                sigman = sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)),
                cin.l = round(pmax(N*psiin - 1.96*sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)), N)),
                cin.u = round(N*psiin + 1.96 *sqrt(N^2*sigma^2 + N*psiin*(psiin - 1)))) %>% as.data.frame()
    result = subset(result, select = -c("psiin"))
    object = list(result = result, N = N)
    class(object) = "psinhat"
    return(object)
  }else{

    #renaming the columns of List_matrix for ease of use
    colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

    if(nfolds > 1 & nfolds > N/50) {
      nfolds = pmax(floor(N/50), 1)
      cat("nfolds is reduced to ", nfolds, " to have sufficient test data.\n")
    }
    psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3*length(funcname))
    rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(k) {
      sapply((k + 1):K, function(s) {
        return(paste0(k, ",", s))
      })}))
    colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
    var_summary = psiinv_summary

    ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname) + 1)
    colnames(ifvals) = c("listpair", funcname)
    ifvals[,"listpair"] = rep(rownames(psiinv_summary), each = N)

    nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname) + 1)
    colnames(nuis) = c("listpair", paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.'))
    nuis[,"listpair"] = ifvals[,"listpair"]
    nuistmle = nuis

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
        psiinvmat = matrix(numeric(0), nrow = nfolds, ncol = 3*length(funcname))
        colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
        varmat = psiinvmat

        ifvalsfold = matrix(numeric(0), nrow = N, ncol = length(funcname))
        colnames(ifvalsfold) = funcname

        nuisfold = matrix(numeric(0), nrow = N, ncol = 3*length(funcname))
        colnames(nuisfold) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
        nuistmlefold = nuisfold

        for(folds in 1:nfolds){#print(folds)

          if(nfolds == 1){
            List1 = List_matrix
            List2 = List1
            sbset = 1:N
          }else{
            sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
            sbset = sbset[sbset <= N]
            List1 = List_matrix[permutset[-sbset],]
            List2 = List_matrix[permutset[sbset],]
          }

          yi = List2[,paste("L", i, sep = '')]
          yj = List2[,paste("L", j, sep = '')]

          overlapij = mean(List1[,i]*List1[,j])
          if(overlapij < eps) {
            warning(cat("Overlap between the lists ", i, " and ", j, " is less than ", eps, '.\n', sep = ''))
          }
          for (func in funcname){

            #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
            qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, i, j, eps = eps, ...), silent = TRUE)

            if ("try-error" %in% class(qhat)) {
              next
            }

            q12 = qhat$q12
            q1 = pmin(pmax(q12, qhat$q1), 1)
            q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))

            nuisfold[sbset, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

            gammainvhat = q1*q2/q12
            psiinvhat = mean(gammainvhat, na.rm = TRUE)

            phihat = gammainvhat*(yj/q2 + yi/q1 - yi*yj/q12) - psiinvhat
            ifvalsfold[sbset, func] = phihat

            Qnphihat = mean(phihat, na.rm = TRUE)

            psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

            psiinvmat[folds, paste(func, c("PI", "DR"), sep = '.')] = c(psiinvhat, psiinvhat.dr)

            sigmasq = var(phihat, na.rm = TRUE)
            varmat[folds, paste(func, c("PI", "DR"), sep = '.')] = sigmasq/N

            datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12))
            colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")

            if(TMLE) {
              tmle = tmle(datmat = datmat, eps = eps, K = K, ...)
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

              nuistmlefold[sbset, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

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

        ifvals[ifvals[,"listpair"] == paste0(i, ",", j), colnames(ifvals) != "listpair"] = ifvalsfold
        nuis[nuis[, "listpair"] == paste0(i, ",", j), colnames(nuis) != "listpair"] = nuisfold
        nuistmle[nuistmle[, "listpair"] == paste0(i, ",", j), colnames(nuistmle) != "listpair"] = nuistmlefold
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
    object = list(result = result, N = N, ifvals = as.data.frame(ifvals), nuis = as.data.frame(nuis), nuistmle = as.data.frame(nuistmle))
    class(object) = "psinhat"
    #print.psinhat(result)
    return(object)
  }
}
#' @export
print.psinhat <- function(obj){#} = "psinhat"){
  print(obj$result)
  invisible(obj)
}
