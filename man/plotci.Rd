% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotci.R
\name{plotci}
\alias{plotci}
\title{Estimate total population size and capture probability using user provided set of models.}
\usage{
plotci(psinhat, psinhatcond, show.plot = TRUE, tsize = 12)
}
\arguments{
\item{psinhat}{An object of type psinhat returned value.}

\item{psinhatcond}{An object of type psinhat returned value.}

\item{show.plot}{A logical value indicating whether it will show plots.}

\item{tsize}{The text size for the plots.}
}
\value{
A list of containing the following components:
\item{result}{  A dataframe of the values in \code{psinhat} which can be passed to ggplot.}
\item{sigma2}{  A dataframe of the values in \code{psinhatcond} which can be passed to ggplot.}
\item{g1}{  A ggplot object with population size estimates and the 95\% confidence interval of the population size \code{n}.}
\item{g2}{  A ggplot object with population size estimates and the 95\% confidence interval of the population size \code{n} conditional on \code{condvar}.}
}
\description{
Estimate total population size and capture probability using user provided set of models.
}
\examples{
n = 10000
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))

p = psinhat(List_matrix = datacrc, funcname = c("logit", "gam", "sl"))
plotci(psinhat = p)

ss = sample(c('a','b','c','d','e','f'), nrow(datacrc), replace = TRUE, prob = (1:6)/sum(1:6))

datacrc1 = data.frame(datacrc, ss)
p = psinhatcond(List_matrix = datacrc1, condvar = 'ss')
plotci(psinhatcond = p)
}
