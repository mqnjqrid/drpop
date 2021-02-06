#' Estimate total population size and capture probability using user provided set of models.
#'
#' @param List_matrix The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists in the data. typically the first \code{K} rows of List_matrix.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param twolist The logical value of whether targeted maximum likelihood algorithm fits only two modes when K = 2.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param iter An integer denoting the maximum number of iterations allowed for targeted maximum likelihood method.
#' @param sl.lib algorithm library for SuperLearner. Default library includes "gam", "glm", "glmnet", "glm.interaction", "ranger".
#' @return A list of estimates containing the following components:
#' \item{psi}{  A dataframe of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The rows represent the list pair which is assumed to be independent conditioned on the covariates.
#' The columns represent the model and method combinations (PI = plug-in, BC = bias-corrected, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
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
psinhat <- function(List_matrix, K = 2, filterrows = TRUE, funcname = c("logit"), nfolds = 5, twolist = FALSE, eps = 0.005, iter = 50, sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet")){

  l = ncol(List_matrix) - K
  n = nrow(List_matrix)

  stopifnot(!is.null(dim(List_matrix)))

  stopifnot(informat(List_matrix = List_matrix, K = K))

  List_matrix = na.omit(List_matrix)

  if(filterrows){
    #removing all rows with only 0's
    List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  }

  List_matrix = as.data.frame(List_matrix)
  #N = number of observed or captured units
  N = nrow(List_matrix)

  stopifnot(((l == 0)&(nrow(List_matrix)>50)) | (nrow(List_matrix) > 0))

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

    psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 1)
    rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(k) {
      sapply((k + 1):K, function(s) {
        return(paste(k, ", ", s, sep = ''))
      })}))
    colnames(psiinv_summary) =  c("PI")
    var_summary = psiinv_summary
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
        psiinv_summary[paste(i, ", ", j, sep = ''),] = q1*q2/q12
        var_summary[paste(i, ", ", j, sep = ''),] = q1*q2*(q1*q2 - q12)*(1 - q12)/q12^3/N
        ifvals = NULL
      }
    }
    return(list(psi = 1/psiinv_summary, sigma2 = N*var_summary, n = N*psiinv_summary,
                varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1), N = N
    ))
  }else{

    #converting factor columns to numeric
    List_matrix = reformat(List_matrix)
    #renaming the columns of List_matrix for ease of use
    colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

    if(nfolds > 1 & nfolds > N/50) {
      nfolds = pmax(floor(N/50), 1)
      cat("nfolds is reduced to ", nfolds, " to have sufficient training data.\n")
    }
    psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3*length(funcname))
    rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(k) {
      sapply((k + 1):K, function(s) {
        return(paste(k, ", ", s, sep = ''))
      })}))
    colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "BC", "TMLE"), sep = '.')
    var_summary = psiinv_summary

    ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname))
    colnames(ifvals) = funcname
    rownames(ifvals) = rep(rownames(psiinv_summary), each = N)

    nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname))
    colnames(nuis) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
    rownames(nuis) = rownames(ifvals)
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
        psiinvmat = matrix(NA, nrow = nfolds, ncol = 3*length(funcname))
        colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "BC", "TMLE"), sep = '.')
        varmat = psiinvmat

        ifvalsfold = matrix(NA, nrow = N, ncol = length(funcname))
        colnames(ifvalsfold) = funcname

        nuisfold = matrix(NA, nrow = N, ncol = 3*length(funcname))
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

          if(mean(List1[,i]*List1[,j]) > eps) {

            for (func in funcname){

              colsubset = stringr::str_subset(colnames(psiinv_summary), func)
              qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, i, j, eps, sl.lib = sl.lib), silent = TRUE)

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

              psiinvhatq = max(psiinvhat + Qnphihat, 1)

              psiinvmat[folds, colsubset][1:2] = c(psiinvhat, psiinvhatq)

              sigmasq = var(phihat, na.rm = TRUE)
              varmat[folds, colsubset][1:2] = sigmasq/N

              datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12))
              datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
              colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")

              tmle = tmle(datmat = datmat, iter = iter, eps = eps, eps_stop = 0.00001, twolist = twolist, K = K)

              if(tmle$error){
                warning("TMLE did not run or converge.")
                psiinvmat[folds,colsubset][3] = NA
                varmat[folds,colsubset][3] = NA
              }else{
                datmat = tmle$datmat
                q12 = pmax(datmat$q12, eps)
                q1 = pmin(datmat$q12 + datmat$q10, 1)
                q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

                nuistmlefold[sbset, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

                gammainvhat = q1*q2/q12
                psiinvhat = mean(gammainvhat, na.rm = TRUE)

                phihat = gammainvhat*(yi/q1 + yj/q2 - yi*yj/q12) - psiinvhat

                Qnphihat = mean(phihat, na.rm = TRUE)

                psiinvmat[folds,colsubset][3] = psiinvhat
                sigmasq = var(phihat, na.rm = TRUE)
                varmat[folds,colsubset][3] = sigmasq/N
              }
            }
          }else{
            message(cat("Overlap between the lists", i, "and", j, "is less than", eps))
            psiinvmat[folds,] = NA
            varmat[folds,] = NA
          }
        }

        psiinv_summary[paste(i, ", ", j, sep = ''),] = colMeans(psiinvmat, na.rm = TRUE)
        var_summary[paste(i, ", ", j, sep = ''),] = colMeans(varmat, na.rm = TRUE)

        ifvals[rownames(ifvals) == paste(i, ", ", j, sep = ''),] = ifvalsfold
        nuis[rownames(nuis) == paste(i, ", ", j, sep = ''),] = nuisfold
        nuistmle[rownames(nuistmle) == paste(i, ", ", j, sep = ''),] = nuistmlefold
      }
    }

    return(list(psi = 1/psiinv_summary, sigma2 = N*var_summary, n = N*psiinv_summary,
                varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1), N = N,
                ifvals = ifvals, nuis = nuis, nuistmle = nuistmle,
                cin.l = pmax(N*psiinv_summary - 1.96*sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)), N),
                cin.u = N*psiinv_summary + 1.96 *sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1))
    ))
  }
}
