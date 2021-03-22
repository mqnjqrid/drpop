#' Estimate the nuisance functions using user provided set of models.
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
#' @param Nmin The cutoff for minimum sample size to perform doubly robust estimation. Otherwise, Petersen estimator is returned.
#' @param num_cores The number of cores to be used for paralellization in Super Learner.
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
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2,1), nrow = nrow(data))
#'
#' data = cbind(data, x)
#' qhat_estimate = qhateval(List_matrix = data, funcname = c("logit", "sl"), nfolds = 2, twolist = FALSE, eps = 0.005)
#' @export
qhateval <- function(List_matrix, K = 2, filterrows = TRUE, funcname = c("logit"), nfolds = 5, twolist = FALSE, eps = 0.005, iter = 50,
                    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.glmnet"), Nmin = 500, num_cores = NA){
  
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
    print("No covariates or not sufficient data.")
    return(0)
  }else{
    
    #converting factor columns to numeric
    List_matrix = reformat(List_matrix)
    #renaming the columns of List_matrix for ease of use
    colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))
    
    if(nfolds > 1 & nfolds > N/50) {
      nfolds = pmax(floor(N/50), 1)
      cat("nfolds is reduced to ", nfolds, " to have sufficient training data.\n")
    }
    
    q1mat = matrix(NA, nrow = N, ncol = length(funcname)*nfolds)
    colnames(q1mat) = paste(rep(funcname, nfolds), rep(1:nfolds, each = length(funcname)), sep = '.')
    q1mat = as.data.frame(q1mat)
    q2mat = q1mat
    q12mat = q1mat
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

        for(folds in 1:nfolds){#print(folds)
          
          sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
          sbset = sbset[sbset <= N]
          List1 = List_matrix[permutset[-sbset],]
          List2 = List_matrix
          
          yi = List2[,paste("L", i, sep = '')]
          yj = List2[,paste("L", j, sep = '')]
          
          if(mean(List1[,i]*List1[,j]) > eps) {
            
            for (func in funcname){
              
              #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
              qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, i, j, eps, sl.lib = sl.lib, num_cores = num_cores), silent = TRUE)
              
              if ("try-error" %in% class(qhat)) {
                next
              }
              
              q12mat[,paste0(func, '.', folds)] = qhat$q12
              q1mat[,paste0(func, '.', folds)] = qhat$q1
              q2mat[,paste0(func, '.', folds)] = qhat$q2
            }
          }else{
            message(cat("Overlap between the lists", i, "and", j, "is less than", eps))
          }
        }
      }
    }
    
    return(list(q1mat = q1mat,
                q2mat = q2mat,
                q12mat = q12mat
    ))
  }
}

