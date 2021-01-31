#' Estimate total population size and capture probability using user provided set of models conditioned on an attribute.
#'
#' @param List_matrix The data frame in capture-recapture format for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param K The number of lists in the data. typically the first \code{K} rows of List_matrix.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param condvar The covariate for which conditional estimates are required.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param twolist The logical value of whether targeted maximum likelihood algorithm fits only two modes when K = 2.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param iter An integer denoting the maximum number of iterations allowed for targeted maximum likelihood method.
#' @param sl.lib algorithm library for SuperLearner. Default library includes "gam", "glm", "glmnet", "glm.interaction", "ranger".
#' @return A list of estimates containing the following components:
# \item{psiinvmat}{ A dataframe of estimated psi inverse for the folds, list pair, model and method combination.
#      The \code{listpair} column represents the list pair which is assumed to be independent conditioned on the covariates.
#      The columns represent the model and method combinations (PI = plug-in, BC = bias-corrected, TMLE = targeted maximum likelihood estimate).}
# \item{varmat}{ A dataframe of estimated sigma^2 in the same format as \code{psimat}.}
#' \item{psi}{  A dataframe of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The \code{listpair} column represents the list pair which is assumed to be independent conditioned on the covariates.
#' The \code{condvar} column stores the value of the original \code{condvar} attribute in \code{List_matrix} which is conditioned upon.
#' The columns represent the model and method combinations (PI = plug-in, BC = bias-corrected, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
#' \item{sigma2}{  A dataframe of the efficiency bound \code{sigma^2} in the same format as \code{psi}.}
#' \item{n}{  A dataframe of the estimated population size n in the same format as \code{psi}.}
#' \item{varn}{  A dataframe of the variance for population size estimate in the same format as \code{psi}.}
#' \item{N}{  The number of data points used in the estimation for each value of \code{condvar}.}
#' \item{cin.l}{  The daraframe of estimated lower bound of a 95% confidence interval of \code{n}.}
#' \item{cin.u}{  The dataframe of estimated upper bound of a 95% confidence interval of \code{n}.}
#'
#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2,1), nrow = nrow(data))
#' ss = sample(1:6, nrow(data), replace = TRUE)
#'
#' data = cbind(data, x, ss)
#' psin_estimate = psinhatcond(List_matrix = data, funcname = c("logit", "sl"), condvar = 'ss', nfolds = 2, twolist = FALSE, eps = 0.005)
#' #this returns the plug-in, the bias-corrected and the tmle estimate for the two models conditioned on column ss
#' @export
psinhatcond <- function(List_matrix, K = 2, filterrows = FALSE, funcname = c("logit"), condvar, nfolds = 5, twolist = FALSE, eps = 0.005, iter = 50, sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet")){

  l = ncol(List_matrix) - K
  n = nrow(List_matrix)

  stopifnot(!is.null(dim(List_matrix)))

  stopifnot(!missing(condvar))
  stopifnot(is.element(condvar, colnames(List_matrix)))

  List_matrix = as.data.frame(List_matrix)
  #N = number of observed or captured units
  N = nrow(List_matrix)

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
      condvar = which(colnames(List_matrix) == condvar) - K
    }
  }

  condvar_vec = unique(List_matrix[,condvar + K])

  psi = numeric(0)
  sigma2 = numeric(0)
  n = numeric(0)
  varn = numeric(0)
  N = numeric(0)
  cin.l = numeric(0)
  cin.u = numeric(0)

  for(cvar in condvar_vec){

    List_matrixsub = List_matrix[List_matrix[,K + condvar] == cvar, -c(K + condvar)]
    est = try(psinhat(List_matrix = List_matrixsub, K = K, filterrows = filterrows, funcname = funcname, nfolds = 2, twolist = twolist, eps = eps, iter = iter, sl.lib = sl.lib), silent = TRUE)

    if("try-error" %in% class(est)){
      next
    }

    psi = rbind(psi, data.frame(listpair = rownames(est$psi), est$psi, condvar = cvar), make.row.names = FALSE)
    sigma2 = rbind(sigma2, data.frame(listpair = rownames(est$psi), est$sigma2, condvar = cvar), make.row.names = FALSE)
    n = rbind(n, data.frame(listpair = rownames(est$psi), est$n, condvar = cvar), make.row.names = FALSE)
    varn = rbind(varn, data.frame(listpair = rownames(est$psi), est$varn, condvar = cvar), make.row.names = FALSE)
    N = rbind(N, data.frame(N = est$N, condvar = cvar))
    cin.l = rbind(cin.l, data.frame(listpair = rownames(est$psi), est$cin.l, condvar = cvar), make.row.names = FALSE)
    cin.u = rbind(cin.u, data.frame(listpair = rownames(est$psi), est$cin.u, condvar = cvar), make.row.names = FALSE)
  }
  return(list(psi = psi, sigma2 = sigma2, n = n, varn = varn, N = N, cin.l = cin.l, cin.u = cin.u))
}
