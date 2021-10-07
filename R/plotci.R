#' Plot estimated confidence interval of total population size.
#'
#' @param object An object of class \code{popsize} or \code{popsize_cond}.
#' @param tsize The text size for the plots.
#' @param ... Any extra arguments passed into the function.
#' @return A ggplot object \code{fig} with population size estimates and the 95% confidence intervals.
#' @examples
#'
#' data = simuldata(n = 10000, l = 1)$data_xstar
#'
#' p = popsize(data = data, funcname = c("logit", "gam"))
#' plotci(p)
#'
#' data = simuldata(n = 10000, l = 1, categorical = TRUE)$data_xstar
#
#' p = popsize_cond(data = data, condvar = 'catcov')
#' plotci(p)
#' @import ggplot2
#' @export
plotci <- function(object, tsize = 12,...){
  stopifnot(!missing(object))
  require(ggplot2, quietly = TRUE)
  fig = NA
  if(class(object) == "popsize"){

    result <- object$result

    fig <- ggplot(result, aes(x = model, color = method)) +
      #geom_line(aes(y = n, linetype = method)) +
      geom_point(aes(y = n), position=position_dodge(0.35)) +
      geom_errorbar(aes(ymin = cin.l, ymax = cin.u), width=.2, position=position_dodge(0.35)) +
      facet_wrap(~listpair, labeller = label_both) +
      scale_color_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9")) +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = tsize))

  }else if(class(object) == "popsize_cond"){

    result <- object$result
    N <- object$N

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
    cat("object is not of class popsize or popsize_cond.\n")
  }
  return(fig)
}
