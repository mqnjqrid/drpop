#' Estimate total population size and capture probability using user provided set of models conditioned on an attribute.
#'
#' @param List_matrix The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists in the data. typically the first \code{K} rows of List_matrix.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param condvar The covariate for which conditional estimates are required.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param TMLE The logical value to indicate whether TMLE has to be computed.
#' @param PLUGIN The logical value to indicate whether the plug-in estimates is returned.
#' @return A list of estimates containing the following components for each list-pair, model and method (PI = plug-in, DR = doubly-robust, TMLE = targeted maximum likelihood estimate):
#' \item{result}{  A dataframe of the below estimated quantities.
#' \itemize{
#' \item{psi}{  The estimated capture probability.}
#' \item{sigma}{  The efficiency bound.}
#' \item{n}{  The estimated population size n.}
#' \item{sdn}{  The estimated standard deviation of the population size.}
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
#' ss = sample(1:6, nrow(data), replace = TRUE)
#'
#' data = cbind(data, x, ss)
#' psin_estimate = popsize_cond(List_matrix = data, funcname = c("logit", "sl"), condvar = 'ss', nfolds = 2, twolist = FALSE, eps = 0.005)
#' #this returns the plug-in, the bias-corrected and the tmle estimate for the two models conditioned on column ss
#' @export
popsize_cond <- function(List_matrix, K = 2, filterrows = FALSE, funcname = c("rangerlogit"), condvar, nfolds = 2, eps = 0.005, TMLE = TRUE, PLUGIN = TRUE, ...){

  l = ncol(List_matrix) - K
  n = nrow(List_matrix)

  stopifnot(!is.null(dim(List_matrix)))

  stopifnot(!missing(condvar))
  stopifnot(is.element(condvar, colnames(List_matrix)))

  List_matrix = as.data.frame(List_matrix)

  conforminglists = apply(List_matrix[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    stop("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }

  if(sum(conforminglists) < K){
    Message(cat("Lists ", which(conforminglists == FALSE), " are not in the required format."))
  }

  if(!missing(condvar)){
    if(is.character(condvar)){
      stopifnot(condvar %in% colnames(List_matrix))
      condvar = which(colnames(List_matrix) == condvar) - K
    }
  }

  condvar_vec = unique(List_matrix[, condvar + K])

  object = NULL

  for(cvar in condvar_vec){

    List_matrixsub = List_matrix[List_matrix[,K + condvar] == cvar, -c(K + condvar)]
    est = try(popsize(List_matrix = List_matrixsub, K = K, filterrows = filterrows, funcname = funcname, nfolds = nfolds, ...), silent = TRUE)

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
print.popsize_cond <- function(obj){
  obj$result$psi = round(obj$result$psi, 3)
  obj$result$sigma = round(obj$result$sigma, 3)
  obj$result$sdn = round(obj$result$sdn, 3)
  print(obj$result)
  invisible(obj)
}
