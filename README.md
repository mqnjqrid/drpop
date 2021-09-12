
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drpop

<!-- badges: start -->
<!-- badges: end -->

The goal of drpop is to provide users doubly robust and efficient
estimates of population size and the condifence intervals for a
capture-recapture problem.

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

## Examples

This is a basic example which shows you how to solve a common problem:

``` r
library(drpop)
#> Registered S3 methods overwritten by 'tibble':
#>   method     from  
#>   format.tbl pillar
#>   print.tbl  pillar

n = 1000
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))[y1+y2 > 0, ]

estim <- popsize(List_matrix = datacrc, func = c("gam"), nfolds = 2, K = 2)
#> Warning: package 'dplyr' was built under R version 4.0.5
#> Loading required package: tidyr
#> Loaded gam 1.20
#> Warning in popsize_base(List_matrix, K = K, j0 = j, k0 = k, filterrows =
#> filterrows, : Plug-in variance is not well-defined. Returning variance evaluated
#> using DR estimator formula
# The population size estimates are 'n' and the standard deviations are 'sigman'
print(estim)
#>   listpair model method   psi sigma    n sigman cin.l cin.u
#> 1      1,2   gam     DR 0.797 1.282  973 39.008   896  1049
#> 2      1,2   gam     PI 0.793 1.282  978 39.112   901  1055
#> 3      1,2   gam   TMLE 0.760 2.082 1020 60.687   901  1139
## basic example code
```

The following shows the confidence interval of estimates for a toy data.
Real population size is 3000.

``` r
library(drpop)

n = 3000
x = matrix(rnorm(n*3, 2,1), nrow = n)

expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))[y1+y2>0,]

estim <- popsize(List_matrix = datacrc, func = c("gam", "rangerlogit"), nfolds = 2, eps = 0.01)
#> Warning in popsize_base(List_matrix, K = K, j0 = j, k0 = k, filterrows =
#> filterrows, : Plug-in variance is not well-defined. Returning variance evaluated
#> using DR estimator formula
print(estim)
#>   listpair       model method   psi sigma    n sigman cin.l cin.u
#> 1      1,2         gam     DR 0.793 0.989 3034 56.096  2924  3144
#> 2      1,2         gam     PI 0.801 0.989 3004 55.693  2895  3113
#> 3      1,2         gam   TMLE 0.801 0.957 3001 54.295  2894  3107
#> 4      1,2 rangerlogit     DR 0.800 0.742 3008 45.574  2918  3097
#> 5      1,2 rangerlogit     PI 0.824 0.742 2919 44.126  2832  3005
#> 6      1,2 rangerlogit   TMLE 0.810 0.915 2970 52.048  2868  3072
plotci(estim)
#> Warning: package 'ggplot2' was built under R version 4.0.4
```

<img src="man/figures/README-plot1-1.png" width="100%" />

``` r
#result = melt(as.data.frame(estim), variable.name = "estimator", value.name = "population_size")
#ggplot(result, aes(x = population_size - n, fill = estimator, color = estimator)) +
#  geom_density(alpha = 0.4) +
#  xlab("Bias on n")
```

The following shows confidence interval of estimates for toy data with a
categorical covariate. Real population size is 10000.

``` r
library(drpop)
n = 10000
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
catcov = sample(c('m','f'), n, replace = TRUE, prob = c(0.45, 0.55))

y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = sapply(1:n, function(i) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.3*(catcov[i] == 'm') + 0.3*x[i,]), expit(-0.6 + 0.3*(catcov[i] == 'm') + 0.3*x[i,])))})
datacrc = cbind.data.frame(y1, y2, exp(x/2), catcov)[y1+y2>0,]

result = popsize_cond(List_matrix = datacrc, condvar = 'catcov')
#> Warning in popsize_base(List_matrix = List_matrixsub, K = K, filterrows =
#> filterrows, : Plug-in variance is not well-defined. Returning variance evaluated
#> using DR estimator formula

#> Warning in popsize_base(List_matrix = List_matrixsub, K = K, filterrows =
#> filterrows, : Plug-in variance is not well-defined. Returning variance evaluated
#> using DR estimator formula
fig = plotci(result)
fig + geom_hline(yintercept = table(catcov), color = "brown", linetype = "dashed")
```

<img src="man/figures/README-plot2-1.png" width="100%" />
