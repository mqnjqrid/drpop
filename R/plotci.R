#' Plot estimated confidence interval of total population size.
#'
#' @param object An object of class \code{psinhat} or \code{psinhatcond}.
#' @param show.plot A logical value indicating whether it will show plots.
#' @param tsize The text size for the plots.
#' @return A list of containing the following components:
#' \item{result}{  A dataframe of the values in \code{object} which can be passed to ggplot.}
#' \item{fig}{  A ggplot object with population size estimates and the 95% confidence interval of the population size \code{n}.}
#' @examples
#' n = 10000
#' x = matrix(rnorm(n*3, 2, 1), nrow = n)
#' expit = function(xi) {
#'   exp(sum(xi))/(1 + exp(sum(xi)))
#' }
#' y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
#' y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
#' datacrc = cbind(y1, y2, exp(x/2))
#'
#' p = psinhat(List_matrix = datacrc, funcname = c("logit", "gam"))
#' plotci(p)
#'
#' ss = sample(c('a','b','c','d','e','f'), nrow(datacrc), replace = TRUE, prob = (1:6)/sum(1:6))
#'
#' datacrc1 = data.frame(datacrc, ss)
#
#' p = psinhatcond(List_matrix = datacrc1, condvar = 'ss')
#' plotci(p)
#' @export
plotci <- function(object, show.plot = TRUE, tsize = 12, ...){
  stopifnot(!missing(object))
  require(ggplot2, quietly = TRUE)
  require(reshape2, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  result = NA
  fig = NA
  if(class(object) == "psinhat"){
    psi <- reshape2::melt(object$psi, value.name = "psi") %>% separate(Var2, c("model", "method"))
    sigma2 <- reshape2::melt(object$sigma2, value.name = "sigma2") %>% separate(Var2, c("model", "method"))
    n <- reshape2::melt(object$n, value.name = "n") %>% separate(Var2, c("model", "method"))
    varn <- reshape2::melt(object$varn, value.name = "varn") %>% separate(Var2, c("model", "method"))
    N <- object$N
    cin.l <- reshape2::melt(object$cin.l, value.name = "cin.l") %>% separate(Var2, c("model", "method"))
    cin.u <- reshape2::melt(object$cin.u, value.name = "cin.u") %>% separate(Var2, c("model", "method"))

    result<- merge(psi, sigma2, by = c("Var1", "model", "method")) %>%
      merge(n, by = c("Var1", "model", "method")) %>%
      merge(varn, by = c("Var1", "model", "method")) %>%
      merge(cin.l, by = c("Var1", "model", "method")) %>%
      merge(cin.u, by = c("Var1", "model", "method")) %>% dplyr::rename(listpair = Var1)
    result$method <- factor(result$method, levels = c("PI", "DR", "TMLE"))

    result <- na.omit(result)

    fig <- ggplot(result, aes(x = model, color = method)) +
      #geom_line(aes(y = n, linetype = method)) +
      geom_point(aes(y = n), position=position_dodge(0.35)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u), width=.2, position=position_dodge(0.35)) +
      facet_wrap(~listpair, labeller = label_both) +
      scale_color_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9")) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = tsize))
  }else if(class(object) == "psinhatcond"){
    psi <- reshape2::melt(object$psi, id.vars = c("listpair", "condvar"), value.name = "psi") %>% separate(variable, c("model", "method"))
    sigma2 <- reshape2::melt(object$sigma2, id.vars = c("listpair", "condvar"), value.name = "sigma2") %>% separate(variable, c("model", "method"))
    n <- reshape2::melt(object$n, id.vars = c("listpair", "condvar"), value.name = "n") %>% separate(variable, c("model", "method"))
    varn <- reshape2::melt(object$varn, id.vars = c("listpair", "condvar"), value.name = "varn") %>% separate(variable, c("model", "method"))
    N <- object$N
    cin.l <- reshape2::melt(object$cin.l, id.vars = c("listpair", "condvar"), value.name = "cin.l") %>% separate(variable, c("model", "method"))
    cin.u <- reshape2::melt(object$cin.u, id.vars = c("listpair", "condvar"), value.name = "cin.u") %>% separate(variable, c("model", "method"))

    result <- merge(psi, sigma2, by = c("listpair", "condvar", "model", "method")) %>%
             merge(n, by = c("listpair", "condvar", "model", "method")) %>%
             merge(varn, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.l, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.u, by = c("listpair", "condvar", "model", "method")) %>%
             merge(N, by = "condvar")
    result$method <- factor(result$method, levels = c("PI", "DR", "TMLE"))

    result <- na.omit(result)

    fig <- ggplot(result, aes(x = condvar, color = method)) +
      #geom_line(aes(y = n, linetype = method)) +
      geom_point(aes(y = n), position=position_dodge(0.35)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u), width=.2, position=position_dodge(0.35)) +
      facet_grid(listpair~model, labeller = label_both) +
      scale_x_discrete(name = "conditional variable (number of observations)", breaks = c(N$condvar), labels = paste(N$condvar, " (", N$N, ')', sep = '')) +
      scale_color_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9")) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = tsize))
  }else{
    cat("object not of class psinhat or psinhatcond\n")
  }
  if(show.plot){
    plot(fig)
  }
  return(invisible(list(result = result, fig = fig)))
}
