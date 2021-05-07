#' A function to reorder the columns of a data table/matrix/data frame and to change factor variables to numeric.
#'
#' @param List_matrix The data table/matrix/data frame which is to be checked.
#' @param capturelists The vector of column names or locations for the capture history list columns.
#' @return \code{List_matrix} With reordered columns so that the capture history columns are followed by the rest.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2, 1), nrow = nrow(data))
#'
#' data = cbind(x, data)
#' result<- reformat(List_matrix = data, capturelists = c(4,5))
#' @export
reformat <- function(List_matrix, capturelists){

  require(tidyverse, quietly = TRUE)

  List_matrix = as.data.frame(List_matrix)

  if(!missing(capturelists)){
    stopifnot(length(capturelists) > 1)
    stopifnot(class(capturelists) %in% c("character", "integer"))
    stopifnot(length(capturelists) <= ncol(List_matrix))

    if(class(capturelists) == "integer"){
      List_matrix = List_matrix[,c(capturelists, setdiff(1:ncol(List_matrix), capturelists))]
    }else{
      List_matrix = List_matrix[,c(capturelists, setdiff(colnames(List_matrix), capturelists))]
    }
    for(i in 1:length(capturelists)){
      List_matrix[,i] <- as.numeric(List_matrix[,i])
     }
  }

  return(List_matrix)
}
