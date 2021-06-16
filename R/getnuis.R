#' Estimate the nuisance functions using user provided set of models.
#'
#' @param List_matrix The data frame in capture-recapture format with two lists for which total population is to be estimated.
#'                    The first K columns are the capture history indicators for the K lists. The remaining columns are covariates in numeric format.
#' @param filterrows A logical value denoting whether to remove all rows with only zeroes.
#' @param funcname The vector of estimation function names to obtain the population size.
#' @param nfolds The number of folds to be used for cross fitting.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param Nmin The cutoff for minimum sample size to perform doubly robust estimation. Otherwise, Petersen estimator is returned.
#' @return A list of estimates containing the following components:
#' \item{q1mat}{  A matrix of the estimated capture probabilities for list1 for each model.}
#' \item{q2mat}{  A matrix of the estimated capture probabilities for list 2.}
#' \item{q12mat}{  A matrix of the estimated capture probabilities for list1 and 2 simultaneously.}
#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @examples
#' data = simuldata(1000, 1)$List_matrix
#' qhat_estimate = getnuis(List_matrix = data, funcname = c("logit", "gam"), nfolds = 2, eps = 0.005)
#' @export
getnuis <- function(List_matrix, i = 1, j = 2, K = 2, filterrows = FALSE, funcname = c("rangerlogit"), nfolds = 5, eps = 0.005, Nmin = 500, ...){

  l = ncol(List_matrix) - K
  n = nrow(List_matrix)

  if(missing(i)){
    i = 1
  }
  if(missing(j)){
    j = 2
  }

  stopifnot(!is.null(dim(List_matrix)))

  if(!informat(List_matrix = List_matrix, K = K)){
    List_matrix = reformat(List_matrix = List_matrix, capturelists = 1:K)
  }

  List_matrix = na.omit(List_matrix)

  if(filterrows){
    #removing all rows with only 0's
    List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  }

  List_matrix = as.data.frame(List_matrix)

  if(K > 2){
    List_matrix = List_matrix[,c(i, j, K+1:l)]
    K = 2
  }
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

    q1mat = matrix(numeric(0), nrow = N, ncol = length(funcname))
    colnames(q1mat) = funcname
    q1mat = as.data.frame(q1mat)
    q2mat = q1mat
    q12mat = q1mat
    idfold = rep(1, N)
    permutset = sample(1:N, N, replace = FALSE)

    for(folds in 1:nfolds){#print(folds)

      sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
      sbset = sbset[sbset <= N]

      idfold[permutset[sbset]] = folds

      List1 = List_matrix[permutset[-sbset],]
      List2 = List_matrix[permutset[sbset],]

      yi = List2[,paste("L", i, sep = '')]
      yj = List2[,paste("L", j, sep = '')]

      overlapij = mean(List1[,i]*List1[,j])
      if(overlapij < eps) {
        warning(cat("Overlap between the lists ", i, " and ", j, " is less than ", eps, '.\n', sep = ''))
      }

        for (func in funcname){

          #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
          qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, i, j, eps, ...), silent = TRUE)

          if ("try-error" %in% class(qhat)) {
            next
          }

          q12mat[permutset[sbset], func] = qhat$q12
          q1mat[permutset[sbset], func] = qhat$q1
          q2mat[permutset[sbset], func] = qhat$q2
        }

    }
    result <- list(q1mat = q1mat, q2mat = q2mat,
                q12mat = q12mat, idfold = idfold)
    class(result) = "getnuis"
    return(result)
  }
}
