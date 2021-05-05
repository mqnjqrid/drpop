#' A function to check whether a given data table/matrix/data frame is in the appropriate for drpop.
#'
#' @param List_matrix The data table/matrix/data frame which is to be checked.
#' @param K The number of lists (optional).
#'
#' @return A boolean for whether \code{List_matrix} is in the appropriate format.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2,1), nrow = nrow(data))
#'
#' informat(List_matrix = data)
#' #this returns TRUE
#'
#' data = cbind(data, x)
#' informat(List_matrix = data)
#' #this returns TRUE
#'
#' informat(List_matrix = data, K = 3)
#' #this returns FALSE
#' @export
informat = function(List_matrix, K){

  if(missing(K)){
    K = 2
  }
  stopifnot(K>1)
  stopifnot(K <= ncol(List_matrix))
  result <- prod(apply(List_matrix[,1:K], 2, function(col){return(setequal(col, c(0,1)))})) > 0
  result <- result & setequal(unlist(lapply(List_matrix[,1:K], "class")), "numeric)
  if(result == FALSE){
    print("The matrix is not in the appropriate format. The first K rows are boolean for capture history information.")
  }
  return(result)
}
