#' A function to reorder the columns of a data table/matrix/data frame and to change factor variables to numeric.
#'
#' @param List_matrix The data table/matrix/data frame which is to be checked.
#' @param capturelists The vector of column names or locations for the capture history list columns.
#' @param condvar The column to obtain conditional estimates. This column is not converted to numeric.
#' @return \code{List_matrix} With reordered columns so that the capture history columns are followed by the rest.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2, 1), nrow = nrow(data))
#'
#' data = cbind(x, data)
#' result<- reformat(List_matrix = data, capturelists = c(4,5))
#' @export
reformat <- function(List_matrix, capturelists, condvar){

  require(tidyverse)

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
  }

  factor_cols <- List_matrix %>% select_if(negate(is.numeric)) %>% names()

  if(!missing(condvar)) {
    factor_cols = setdiff(factor_cols, "condvar")
  }

  if(length(factor_cols)) {
    List_matrix = data.frame(List_matrix[,!(names(List_matrix) %in% factor_cols)],
                           model.matrix(formula(paste('~', paste(factor_cols, collapse = '+'))), List_matrix))
  }

  return(List_matrix)
}
