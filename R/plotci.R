#' Estimate total population size and capture probability using user provided set of models.
#'
#' @param psinhat An object of type psinhat returned value.
#' @param psinhatcond An object of type psinhat returned value.
#' @param show.plot A logical value indicating whether it will show plots.
#' @param tsize The text size for the plots.
#' @return A list of containing the following components:
#' \item{result}{  A dataframe of the values in \code{psinhat} which can be passed to ggplot.}
#' \item{sigma2}{  A dataframe of the values in \code{psinhatcond} which can be passed to ggplot.}
#' \item{g1}{  A ggplot object with population size estimates and the 95% confidence interval of the population size \code{n}.}
#' \item{g2}{  A ggplot object with population size estimates and the 95% confidence interval of the population size \code{n} conditional on \code{condvar}.}
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
#' p = psinhat(List_matrix = datacrc, funcname = c("logit", "gam", "sl"))
#' plotci(psinhat = p)
#'
#' ss = sample(c('a','b','c','d','e','f'), nrow(datacrc), replace = TRUE, prob = (1:6)/sum(1:6))
#'
#' datacrc1 = data.frame(datacrc, ss)
#
#' p = psinhatcond(List_matrix = datacrc1, condvar = 'ss')
#' plotci(psinhatcond = p)
#' @export
plotci <- function(psinhat, psinhatcond, show.plot = TRUE, tsize = 12){
  require(ggplot2, quietly = TRUE)
  require(reshape2, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  result = NA
  resultcond = NA
  g1 = NA
  g2 = NA
  if(!missing(psinhat)){
    psi <- reshape2::melt(psinhat$psi, value.name = "psi")%>% separate(Var2, c("model", "method"))
    sigma2 <- reshape2::melt(psinhat$sigma2, value.name = "sigma2")%>% separate(Var2, c("model", "method"))
    n <- reshape2::melt(psinhat$n, value.name = "n")%>% separate(Var2, c("model", "method"))
    varn <- reshape2::melt(psinhat$varn, value.name = "varn")%>% separate(Var2, c("model", "method"))
    N <- psinhat$N
    cin.l <- reshape2::melt(psinhat$cin.l, value.name = "cin.l")%>% separate(Var2, c("model", "method"))
    cin.u <- reshape2::melt(psinhat$cin.u, value.name = "cin.u")%>% separate(Var2, c("model", "method"))

    result<- merge(psi, sigma2, by = c("Var1", "model", "method")) %>%
      merge(n, by = c("Var1", "model", "method")) %>%
      merge(varn, by = c("Var1", "model", "method")) %>%
      merge(cin.l, by = c("Var1", "model", "method")) %>%
      merge(cin.u, by = c("Var1", "model", "method")) %>% dplyr::rename(listpair = Var1)
    result$method <- factor(result$method, levels = c("PI", "DR", "TMLE"))

    result <- na.omit(result)

    g1 <- ggplot(result, aes(x = model, color = method)) +
      #geom_line(aes(y = n, linetype = method)) +
      geom_point(aes(y = n), position=position_dodge(0.25)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u), width=.2, position=position_dodge(0.25)) +
      facet_wrap(~listpair, labeller = label_both) +
      scale_color_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9")) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = tsize))

    if(show.plot){
      g1
    }
  }

  if(!missing(psinhatcond)){
    psi <- reshape2::melt(psinhatcond$psi, id.vars = c("listpair", "condvar"), value.name = "psi")%>% separate(variable, c("model", "method"))
    sigma2 <- reshape2::melt(psinhatcond$sigma2, id.vars = c("listpair", "condvar"), value.name = "sigma2")%>% separate(variable, c("model", "method"))
    n <- reshape2::melt(psinhatcond$n, id.vars = c("listpair", "condvar"), value.name = "n")%>% separate(variable, c("model", "method"))
    varn <- reshape2::melt(psinhatcond$varn, id.vars = c("listpair", "condvar"), value.name = "varn")%>% separate(variable, c("model", "method"))
    N <- psinhatcond$N
    cin.l <- reshape2::melt(psinhatcond$cin.l, id.vars = c("listpair", "condvar"), value.name = "cin.l")%>% separate(variable, c("model", "method"))
    cin.u <- reshape2::melt(psinhatcond$cin.u, id.vars = c("listpair", "condvar"), value.name = "cin.u")%>% separate(variable, c("model", "method"))

    resultcond <- merge(psi, sigma2, by = c("listpair", "condvar", "model", "method")) %>%
             merge(n, by = c("listpair", "condvar", "model", "method")) %>%
             merge(varn, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.l, by = c("listpair", "condvar", "model", "method")) %>%
             merge(cin.u, by = c("listpair", "condvar", "model", "method")) %>%
             merge(N, by = "condvar")
    resultcond$method <- factor(resultcond$method, levels = c("PI", "DR", "TMLE"))

    resultcond <- na.omit(resultcond)

    g2<- ggplot(resultcond, aes(x = condvar, color = method)) +
      #geom_line(aes(y = n, linetype = method)) +
      geom_point(aes(y = n), position=position_dodge(0.25)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u), width=.2, position=position_dodge(0.25)) +
      facet_grid(listpair~model, labeller = label_both) +
      scale_x_discrete(name = "conditional variable (number of observations)", breaks = c(N$condvar), labels = paste(N$condvar, " (", N$N, ')', sep = '')) +
      scale_color_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9")) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = tsize))

    if(show.plot){
      g2
    }
  }
  return(invisible(list(result = result, resultcond = resultcond, g1 = g1, g2 = g2)))
}
