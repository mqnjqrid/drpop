#' A function to reorder the columns of a data table/matrix/data frame and to change factor variables to numeric.
#'
#' @param n The size of the population.
#' @param l The number of continuous covariates.
#' @param categorical A logical value of whether to include a categorical column.
#' @param ep A numeric value to change the list probabilities.
#' @param K The number of lists. Default value is 2. Maximum value is 3.
#' @return A list of estimates containing the following components:
#' \item{List_matrix}{  A dataframe in with \code{K} list capture histories and covariates from a population if true size \code{n} with only observed rows.}
#' \item{List_matrix_xstar}{  A dataframe in with two list capture histories and transformed covariates from a population if true size \code{n} with only observed rows.}
#' \item{psi0}{  The empirical capture probability for the set-up used.}
#' \item{pi1}{  The conditional capture probabilities for list 1.}
#' \item{pi2}{  The conditional capture probabilities for list 2.}
#' \item{pi3}{  The conditional capture probabilities for list 3 when \code{K = 3}.}
#'
#' @examples
#' data = simuldata(n = 1000, l = 2)$List_matrix
#' psi0 = simuldata(n = 10000, l = 2)$psi0
#' @export
simuldata = function(n, l, categorical = FALSE, ep = 0, K = 2){
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  logit = function(x) {
    log(x/(1 - x))
  }

  x = matrix(rnorm(n*l, 3, 1.5), nrow = n, ncol = l)
  colnames(x) = paste0('x', 1:l)

  if(categorical){
    catcov = factor(sample(c('a','b','c'), n, replace = TRUE, prob = c(1/3, 1/3, 1/3)),
                    levels = c('a','b','c'))#categorical covariate
  }else{
    catcov = NA
  }
  pi1 = function(x, catcov = NA) {
    if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
      x = rowSums(x)
    }
    #expit( 1 + ep + catadjust + 0.4*x)
    expit(ep + 0.4*x)
  }

  pi2 = function(x, catcov = NA) {
    if(sum(is.na(catcov))){
      catadjust = 0
    }else{
      catadjust <- (1 - as.numeric(catcov))*0.8
    }
    if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
      x = rowSums(x)
    }
    #expit( -0.5 + ep + catadjust + 0.3*x)
    expit(ep + catadjust + 0.3*x)
  }
  if(K > 2){
    pi3 = function(x, catcov = NA) {
      if(sum(is.na(catcov))){
        catadjust = 0
      }else{
        catadjust <- (1 - as.numeric(catcov))*0.5
      }
      if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
        x = rowSums(x)
      }
      #expit( 1 + ep + catadjust + 0.1*x)
      expit(ep + catadjust + 0.2*x)
    }
  }
  p1 = pi1(x = x, catcov = catcov)
  p2 = pi2(x = x, catcov = catcov)

  y1 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p1[i], p1[i]))})
  y2 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p2[i], p2[i]))})

  if(K > 2){
    p3 = pi3(x, catcov)
    y3 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p3[i], p3[i]))})
  }
  xp = do.call("cbind", lapply(1:ncol(x),
                               function(li){
                                 if(li%%4 == 1){
                                   return(exp(abs(x[,li]-2)/2))
                                 }else if(li%%4 == 2){
                                   return(x[,li]/(1 + exp(x[,li -1])) + 10)
                                 }else if(li%%4 == 3){
                                   return((x[,li]*x[,li-2]/25 + 0.6)^3)
                                 }else{
                                   return((x[,li -2] + x[,li] + 20)^2)
                                 }
                               }))
  List_matrix = cbind.data.frame(y1, y2, x)
  List_matrix_xstar = cbind.data.frame(y1, y2, xp)
  psi0 = 1 -  mean((1 - p1)*(1 - p2))
  if(K > 2){
    List_matrix = cbind.data.frame(y1, y2, y3, x)
    List_matrix_xstar = cbind.data.frame(y1, y2, y3, xp)
    psi0 = 1 -  mean((1 - p1)*(1 - p2)*(1 - p3))
  }
  if(categorical){
    List_matrix$catcov = catcov
    List_matrix_xstar$catcov = catcov
  }

  result = list(List_matrix = List_matrix[rowSums(List_matrix[,1:K])>0,],
                List_matrix_xstar = List_matrix_xstar[rowSums(List_matrix[,1:K])>0,],
                psi0 = psi0, pi1 = pi1, pi2 = pi2)
  if(K > 2){
    result$pi3 = pi3
  }
  return(result)
}
if(FALSE){
simuldata = function(n, l, categorical = FALSE, ep = 0, K = 2){
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  logit = function(x) {
    log(x/(1 - x))
  }

  x = matrix(rnorm(n*l, 2, 1), nrow = n, ncol = l)
  colnames(x) = paste0('x', 1:l)

  if(categorical){
    catcov = factor(sample(c('a','b','c'), n, replace = TRUE, prob = c(1/3, 1/3, 1/3)),
                    levels = c('a','b','c'))#categorical covariate
  }else{
    catcov = NA
  }
  pi1 = function(x, catcov) {
    if(missing(catcov) | is.na(catcov)){
      catadjust = 0
    }else{
      catadjust <- (1 - as.numeric(catcov))*0
    }
    if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
      x = rowSums(x)
    }
    #expit( 1 + ep + catadjust + 0.4*x)
    expit(ep - 0.5 + catadjust + 0.3*x)
  }

  pi2 = function(x, catcov) {
    if(missing(catcov) | is.na(catcov)){
      catadjust = 0
    }else{
      catadjust <- (1 - as.numeric(catcov))*0.8
    }
    if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
      x = rowSums(x)
    }
    #expit( -0.5 + ep + catadjust + 0.3*x)
    expit(ep + catadjust + 0.2*x)
  }
  if(K > 2){
    pi3 = function(x, catcov) {
      if(missing(catcov) | is.na(catcov)){
        catadjust = 0
      }else{
        catadjust <- (1 - as.numeric(catcov))*0.5
      }
      if(sum(!is.null(ncol(x)) & ncol(x) > 1)){
        x = rowSums(x)
      }
      #expit( 1 + ep + catadjust + 0.1*x)
      expit(ep + catadjust - 0.2*x)
    }
  }
  p1 = pi1(x, catcov)
  p2 = pi2(x, catcov)

  y1 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p1[i], p1[i]))})
  y2 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p2[i], p2[i]))})

  if(K > 2){
    p3 = pi3(x, catcov)
    y3 = sapply(1:n, function(i) {sample(c(0, 1), 1, prob = c(1 - p3[i], p3[i]))})
  }
  xp = do.call("cbind", lapply(1:ncol(x),
                               function(li){
                                 if(li%%4 == 1){
                                   return(exp(x[,li]/2))
                                 }else if(li%%4 == 2){
                                   return(x[,li]/(1 + exp(x[,li -1])) + 10)
                                 }else if(li%%4 == 3){
                                   return((x[,li]*x[,li-2]/25 + 0.6)^3)
                                 }else{
                                   return((x[,li -2] + x[,li] + 20)^2)
                                 }
                               }))
  List_matrix = cbind.data.frame(y1, y2, x)
  List_matrix_xstar = cbind.data.frame(y1, y2, xp)
  psi0 = 1 -  mean((1 - p1)*(1 - p2))
  if(K > 2){
    List_matrix$y3 = y3
    List_matrix_xstar$y3 = y3
    psi0 = 1 -  mean((1 - p1)*(1 - p2)*(1 - p3))
  }
  if(categorical){
    List_matrix$catcov = catcov
    List_matrix_xstar$catcov = catcov
  }

  result = list(List_matrix = List_matrix[rowSums(List_matrix[,1:K])>0,],
      List_matrix_xstar = List_matrix_xstar[rowSums(List_matrix[,1:K])>0,],
      psi0 = psi0, pi1 = pi1, pi2 = pi2)
  if(K > 2){
    result$pi3 = pi3
  }
  return(result)
}
}
