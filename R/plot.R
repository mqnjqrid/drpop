#' Estimate total population size and capture probability using user provided set of models.
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
#' @return A list of estimates containing the following components:
#' \item{psi}{  A dataframe of the estimated capture probability for each list pair, model and method combination. In the absence of covariates, the column represents the standard plug-in estimate.
#' The rows represent the list pair which is assumed to be independent conditioned on the covariates.
#' The columns represent the model and method combinations (PI = plug-in, BC = bias-corrected, TMLE = targeted maximum likelihood estimate)indicated in the columns.}
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
#' psin_estimate = psinhat(List_matrix = data)
#' #this returns the basic plug-in estimate since covariates are absent.
#'
#' data = cbind(data, x)
#' psin_estimate = psinhat(List_matrix = data, funcname = c("logit", "sl"), nfolds = 2, twolist = FALSE, eps = 0.005)
#' #this returns the plug-in, the bias-corrected and the tmle estimate for the two models
#' @export
plot <- function(psinhat, psinhatcond){
    require(ggplot2)
    require(reshape2)
    require(tidyr)
  if(!missing(psinhat)){

  }

  if(!missing(psinhatcond)){
    psi <- reshape2::melt(psinhatcond$psi, id.vars = c("listpair", "condvar"), value.name = "psi")%>% separate(variable, c("model", "method"))
    sigma2 <- reshape2::melt(psinhatcond$sigma2, id.vars = c("listpair", "condvar"), value.name = "sigma2")%>% separate(variable, c("model", "method"))
    n <- reshape2::melt(psinhatcond$n, id.vars = c("listpair", "condvar"), value.name = "n")%>% separate(variable, c("model", "method"))
    varn <- reshape2::melt(psinhatcond$varn, id.vars = c("listpair", "condvar"), value.name = "varn")%>% separate(variable, c("model", "method"))
    N <- psinhatcond$N
    cin.l <- reshape2::melt(psinhatcond$cin.l, id.vars = c("listpair", "condvar"), value.name = "cin.l")%>% separate(variable, c("model", "method"))
    cin.u <- reshape2::melt(psinhatcond$cin.u, id.vars = c("listpair", "condvar"), value.name = "cin.u")%>% separate(variable, c("model", "method"))

    result<- merge(psi, sigma2, by = c("listpair", "condvar", "model", "method")) %>%
             merge(n, by = c("listpair", "condvar", "model", "method")) %>%
             merge(varn, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.l, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.u, by = c("listpair", "condvar", "model", "method"))  %>%
             merge(N, by = "condvar")
    ggplot(result, aes(x = condvar, color = method)) +
      geom_line(aes(y = n)) +
      geom_point(aes(y = n)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u)) +
      facet_wrap(~model) +
      theme_bw()
  }
}
