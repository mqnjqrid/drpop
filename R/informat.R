#' A function to check whether a given data table/matrix/data frame is in the appropriate for drpop.
#'
#' @param data The data table/matrix/data frame which is to be checked.
#' @param K The number of lists (optional).
#'
#' @return A boolean for whether \code{data} is in the appropriate format.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2,1), nrow = nrow(data))
#'
#' informat(data = data)
#' #this returns TRUE
#'
#' data = cbind(data, x)
#' informat(data = data)
#' #this returns TRUE
#'
#' informat(data = data, K = 3)
#' #this returns FALSE
#' @export
informat <- function(data, K = 2){

  stopifnot(K>1)
  stopifnot(K <= ncol(data))
  result <- prod(apply(data[,1:K], 2, function(col){return(setequal(col, c(0,1)))})) > 0
  if(!result){
    cat("The matrix is not in the appropriate format.\n")
  }
  result <- result & setequal(unlist(lapply(data[,1:K], "class")), "numeric")
  return(result)
}
