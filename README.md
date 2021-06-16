
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

n = 1000
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c(1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))[y1+y2 > 0, ]

estim <- popsize(List_matrix = datacrc, func = c("gam"), nfolds = 2, K = 2)
#> Warning: package 'tidyverse' was built under R version 4.0.3
#> -- Attaching packages ---------------
#> v ggplot2 3.3.3     v purrr   0.3.4
#> v tibble  3.0.3     v dplyr   1.0.1
#> v tidyr   1.1.1     v stringr 1.4.0
#> v readr   1.3.1     v forcats 0.5.0
#> Warning: package 'ggplot2' was built under R version 4.0.4
#> -- Conflicts ------------------------
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
#> 
#> Attaching package: 'foreach'
#> The following objects are masked from 'package:purrr':
#> 
#>     accumulate, when
#> Loaded gam 1.20
# The population size estimates are 'n' and the standard deviations are 'sigman'
print(estim)
#>   listpair model method   psi sigma   n sigman cin.l cin.u
#> 1      1,2   gam     DR 0.830 0.856 963 27.987   908  1018
#> 2      1,2   gam     PI 0.826 0.856 969 28.127   914  1024
#> 3      1,2   gam   TMLE 0.825 0.764 970 25.945   919  1021
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
print(estim)
#>   listpair       model method   psi sigma    n sigman cin.l cin.u
#> 1      1,2         gam     DR 0.815 0.866 2938 49.622  2841  3035
#> 2      1,2         gam     PI 0.821 0.866 2917 49.318  2821  3014
#> 3      1,2         gam   TMLE 0.834 0.828 2872 47.063  2780  2965
#> 4      1,2 rangerlogit     DR 0.811 0.702 2954 43.242  2870  3039
#> 5      1,2 rangerlogit     PI 0.846 0.702 2833 41.201  2752  2913
#> 6      1,2 rangerlogit   TMLE 0.834 0.729 2872 42.958  2787  2956
plotci(estim)
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
fig = plotci(result)
fig + geom_hline(yintercept = table(catcov), color = "brown", linetype = "dashed")
```

<img src="man/figures/README-plot2-1.png" width="100%" />
