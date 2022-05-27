#' A function to reorder the columns of a data table/matrix/data frame and to change factor variables to numeric.
#'
#' @param data The data table/matrix/data frame which is to be checked.
#' @param capturelists The vector of column names or locations for the capture history list columns.
#' @return \code{data} With reordered columns so that the capture history columns are followed by the rest.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' x = matrix(rnorm(nrow(data)*3, 2, 1), nrow = nrow(data))
#'
#' data = cbind(x, data)
#' result<- reformat(data = data, capturelists = c(4,5))
#' @export
reformat <- function (data, capturelists) 
{
  data = as.data.frame(data)
  if (!missing(capturelists)) {
    stopifnot(length(capturelists) > 1)
    stopifnot(class(capturelists) %in% c("character", "numeric", "integer"))
    stopifnot(length(capturelists) <= ncol(data))
    if (class(capturelists) %in% c("numeric", "integer")) {
      capturelists = round(capturelists)
      data = data[, c(capturelists, setdiff(1:ncol(data), 
                                            capturelists))]
    }
    else {
      data = data[, c(capturelists, setdiff(colnames(data), 
                                            capturelists))]
    }
    for (i in 1:length(capturelists)) {
      data[, i] <- as.numeric(data[, i])
    }
    if(length(capturelists) < ncol(data)) {
      for (i in (ncol(data) - length(capturelists)):ncol(data)) {
        data[,i] <- as.factor(data[,i])
      }
    }
  }
  return(data)
}
