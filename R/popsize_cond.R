#' Estimate total population size and capture probability using user provided set of models conditioned on an attribute.
#'
#' @param data The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists in the data. typically the first \code{K} rows of data.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param condvar The covariate for which conditional estimates are required.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param margin The minimum value the estimates can attain to bound them away from zero.
#' @param sl.lib Algorithm library for [qhat_sl()]. See [SuperLearner::listWrappers()]. Default library includes "gam", "glm", "glmnet", "glm.interaction", "ranger".
#' @param Nmin The cutoff for minimum sample size to perform doubly robust estimation. Otherwise, Petersen estimator is returned.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
#' @param PLUGIN The logical value to indicate whether the plug-in estimates are returned.
#' @param ... Any extra arguments passed into the function. See [qhat_rangerlogit()], [qhat_sl()], [tmle()].
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
#' @seealso \code{\link{popsize}}
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust capture-recapture methods for estimating population size. _arXiv preprint_ *arXiv:2104.14091*.
#' @examples
#' \donttest{
#' data = simuldata(n = 10000, l = 2, categorical = TRUE)$data
#'
#' psin_estimate = popsize_cond(data = data, funcname = c("logit", "gam"),
#'      condvar = 'catcov', PLUGIN = TRUE, TMLE = TRUE)
#' #this returns the plug-in, the bias-corrected and the tmle estimate for the
#' #two models conditioned on column catcov
#' }
#' @export
popsize_cond <- function(data, K = 2, filterrows = FALSE, funcname = c("rangerlogit"), condvar, nfolds = 2, margin = 0.005,
                          sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), TMLE = TRUE, PLUGIN = TRUE, Nmin = 100,...){

  l = ncol(data) - K
  n = nrow(data)

  stopifnot(!is.null(dim(data)))

  stopifnot(!missing(condvar))
  stopifnot(is.element(condvar, colnames(data)))

  data = as.data.frame(data)

  conforminglists = apply(data[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    stop("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }

  if(sum(conforminglists) < K){
    message(cat("Lists ", which(conforminglists == FALSE), " are not in the required format."))
  }

  if(!missing(condvar)){
    if(is.character(condvar)){
      stopifnot(condvar %in% colnames(data))
      condvar = which(colnames(data) == condvar) - K
    }
  }

  condvar_vec = unique(data[, condvar + K])

  object = NULL

  for(cvar in condvar_vec){

    datasub = data[data[,K + condvar] == cvar, -c(K + condvar)]
    est = try(popsize_base(data = datasub, K = K, filterrows = filterrows, funcname = funcname, nfolds = nfolds, margin = margin,
                           sl.lib = sl.lib, Nmin = Nmin, TMLE = TMLE, PLUGIN = PLUGIN, ...), silent = TRUE)

    if("try-error" %in% class(est)){
      next
    }
    if(!('DR' %in% est$result$method)){
      next
    }
    if(is.null(object)){
      object = lapply(est, function(x) cbind.data.frame(x, condvar = cvar))
    }else{
      object = Map("rbind", object, lapply(est, function(x) cbind.data.frame(x, condvar = cvar)))
    }
  }
  if(!is.null(object)){
    colnames(object$N)[1] = 'N'
    class(object) = "popsize_cond"
    return(object)
  }else{
    print("Error in estimation for all subsets.")
    return(0)
  }
}
#' @export
print.popsize_cond <- function(x, ...){
  x$result$psi = round(x$result$psi, 3)
  x$result$sigma = round(x$result$sigma, 3)
  x$result$sigman = round(x$result$sigman, 3)
  print(x$result)
  invisible(x)
}
