
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drpop

<!-- badges: start -->

<!-- badges: end -->

The goal of drpop is to provide users doubly-robust and efficient
estimates of population size and the variances for a capture-recapture
problem.

## Installation

You can install the released version of drpop from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("drpop")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mqnjqrid/drpop")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(drpop)

n = 1000
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))

options(warn = -1)
estim <- psinhat(List_matrix = datacrc, func = c("logit"), nfolds = 2, K = 2)
# The population size estimates are obtained by
estim$n
#>      logit.PI logit.BC logit.TMLE
#> 1, 2 996.2992 981.2471   948.9273
# The corresponding variances are
estim$varn
#>      logit.PI logit.BC logit.TMLE
#> 1, 2 1168.821 1146.431   672.1127
## basic example code
```

The following shows histograms of estimates for toy data. Real
population size is 10000.

``` r
library(drpop)
library(reshape2)
library(ggplot2)

n = 5000
x = matrix(rnorm(n*3, 2,1), nrow = n)

expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))

estim <- do.call("rbind", lapply(1:100, function(iter){
  return(psinhat(List_matrix = datacrc, func = c("logit"), nfolds = 2, eps = 0.01)$n)
}))

result = melt(as.data.frame(estim), variable.name = "estimator", value.name = "population_size")
#> No id variables; using all as measure variables
ggplot(result, aes(x = population_size - n, fill = estimator, color = estimator)) +
  geom_density(alpha = 0.4) +
  xlab("Bias on n")
```

<img src="man/figures/README-plot-1.png" width="100%" />
