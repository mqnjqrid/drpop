#' Returns the targeted maximum likelihood estimates for the nuisance functions
#'
#' @param datmat The data frame containing columns \code{yi}, \code{yj}, \code{yij}, \code{q10}, \code{q02} and \code{q12}.
#' @param iter An integer denoting the maximum number of iterations allowed for targeted maximum likelihood method. Default value is 100.
#' @param eps The minimum value the estimates can attain to bound them away from zero.
#' @param eps_stop The minimum value the estimates can attain to bound them away from zero.
#' @param twolist The logical value of whether targeted maximum likelihood algorithm fits only two modes when K = 2.
#' @param K The number of lists in the original data.
#' @return A list of estimates containing the following components:
#' \item{error}{  An indicator of whether the algorithm ran and converged. Returns FALSE, if it ran correctly and FALSE otherwise.}
#' \item{datmat}{  A data frame returning \code{datmat} with the updated estimates for the nuisance functions \code{q10}, \code{q02} and \code{q12}. This is returned only if \code{error} is FALSE.}

#' @references Gruber, S., & Van der Laan, M. J. (2011). tmle: An R package for targeted maximum likelihood estimation.
#' @examples
#' data = matrix(sample(c(0,1), 2000, replace = TRUE), ncol = 2)
#' xqmat = matrix(runif(nrow(data)*3, 0, 1), nrow = nrow(data))
#' datmat = cbind(data, data[,1]*data[,2], xmat)
#' colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")
#' datmat = as.data.frame(datmat)
#' result = tmle(datmat, eps = 0.005, eps_stop = 0.00001, twolist = TRUE)
#' @export

tmle <- function(datmat, iter = 250, eps = 0.005, eps_stop = 0.005, twolist = FALSE, K = 2,...){

  if(!prod(c("yi", "yj", "yij", "q10", "q02", "q12") %in% colnames(datmat))){
    stop("datmat misses one or more of the following columns: \t (yi, yj, yij, q10, q02, q12).")
    return(list(error = TRUE))
  }

  eps = max(eps, 0.005)
  datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))

 # ifval = (datmat$q10+datmat$q12)*(datmat$q10+datmat$q12)/datmat$q12 *(
#    datmat$yi/(datmat$q10+datmat$q12) + datmat$yj/(datmat$q02+datmat$q12) - datmat$yij/(datmat$q12) - 1)

#  eps_stop = sqrt(mean((ifval)^2))/max(log(nrow(datmat)), 10)/sqrt(nrow(datmat))
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  logit = function(x) {
    log(x/(1 - x))
  }

  epsilon_error = 1 + eps_stop
  cnt = 0

  while (abs(epsilon_error) > eps_stop){
    cnt = cnt + 1
    if (cnt > iter){break}

    ########################### model 1 for q12
    dat1 = cbind(datmat$yij, logit(datmat$q12), (datmat$q10 + datmat$q12)/datmat$q12
                 + (datmat$q02 + datmat$q12)/datmat$q12
                 - (datmat$q10 + datmat$q12)*(datmat$q02 + datmat$q12)/datmat$q12^2 )
    colnames(dat1) = c("yij", "logitq12", "ratio")
    dat1 = as.data.frame(dat1)
    mod1 = try(glm(yij ~ -1 + offset(logitq12) + ratio
                   , family = binomial(link = logit), data = dat1, na.action = na.omit), silent = TRUE)
    if (!("try-error" %in% class(mod1))){
      datmat[,"q12"] = predict(mod1, newdata = dat1, type = "response")

    }
    epsilon_error = abs(mod1$coefficients)
    datmat$q12 = pmax(pmin(datmat$q12, 1), eps)

    ########################### model 2 for q1
    dat1 = cbind(datmat$yi*(1 - datmat$yj), logit(datmat$q10), (datmat$q02 + datmat$q12)/datmat$q12)
    colnames(dat1) = c("yi0", "logitq10", "ratio")
    dat1 = as.data.frame(dat1)
    mod1 = try(glm(yi0 ~ -1 + offset(logitq10) + ratio, family = binomial(link = logit), data = dat1, na.action = na.omit), silent = TRUE)
    if (!("try-error" %in% class(mod1))){
      datmat$q10 = predict(mod1, newdata = dat1, type = "response")
      datmat[,"q10"] = pmin(datmat[,"q10"], 1 - datmat$q12)
    }
    datmat$q10 = pmax(pmin(datmat$q10, 1), eps)
    epsilon_error = max(abs(mod1$coefficients), epsilon_error)

    ########################### model 3 for q2
    if (K > 2 | twolist == FALSE){
      dat1 = cbind(datmat$yj*(1 - datmat$yi), logit(datmat$q02), (datmat$q10 + datmat$q12)/datmat$q12)
      colnames(dat1) = c("y0j", "logitq02", "ratio")
      dat1 = as.data.frame(dat1)
      mod1 = try(glm(y0j ~ -1 + offset(logitq02) + ratio, family = binomial(link = logit), data = dat1, na.action = na.omit), silent = TRUE)
      if (!("try-error" %in% class(mod1))){
        datmat$q02 = predict(mod1, newdata = dat1, type = "response")
        datmat[,"q02"] = pmin(datmat[,"q02"], 1 - datmat$q10 - datmat$q12)
        epsilon_error = max(abs(mod1$coefficients), epsilon_error)
      }
    }else{
      datmat[,"q02"] = pmax(0, 1 - datmat$q10 - datmat$q12)
    }

    datmat$q02 = pmax(pmin(datmat$q02, 1), eps)

    # ifval = (datmat$q10+datmat$q12)*(datmat$q10+datmat$q12)/datmat$q12 *(
    #   datmat$yi/(datmat$q10+datmat$q12) + datmat$yj/(datmat$q02+datmat$q12) - datmat$yij/(datmat$q12) - 1)
    # eplison_error = mean(ifval)
    if(is.null(epsilon_error))
       epsilon_error = 2
    #epsilon_error = max(abs(expit(dat1$logitq12) - datmat$q12),
    #                    abs(expit(dat1$logitq10) - datmat$q10),
    #                    abs(expit(dat1$logitq02) - datmat$q02))
  }
  return(list(error = abs(epsilon_error) > 1, datmat = datmat, iterations = cnt, epsilon_error = epsilon_error))
}
