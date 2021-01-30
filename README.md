
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EffPopSize

<!-- badges: start -->

<!-- badges: end -->

The goal of EffPopSize is to provide users doubly-robust and efficient
estimates of population size and the variances for a capture-recapture
problem.

## Installation

You can install the released version of EffPopSize from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("EffPopSize")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mqnjqrid/EffPopSize")
```

## Example

This is a basic example which shows you how to solve a common problem:

    #>      logit PI logit BC logit TMLE
    #> 1, 2 951.8203 976.6341    966.611
    #>      logit PI logit BC logit TMLE
    #> 1, 2 744.7828 781.1612   793.3336

The following shows histograms of estimates for toy data. Real
population size is 10000.

    #> No id variables; using all as measure variables

<img src="man/figures/README-plot-1.png" width="100%" />
