
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplingNR

<!-- badges: start -->

[![R-CMD-check](https://github.com/jmendelson256/samplingNR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmendelson256/samplingNR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`samplingNR` is an R package that allows for computing optimal
allocations under anticipated nonresponse. The underlying theory is
provided in Mendelson & Elliott (in press).

## Installation

You can install the latest development version of `samplingNR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jmendelson256/samplingNR", build_vignettes = TRUE)
```

## Vignette

You can learn about how to use the package in `vignette("samplingNR")`,
which shows how to replicate two of the tables from our paper in an
application to a post-election survey of military personnel.

## Basic usage

The main allocation function is `opt_nh_nonresp()`, which provides the
exact version of our proposed optimal allocation under anticipated
nonresponse. The current version assumes that the goal is to minimize
the (expected) variance subject to a constraint on the total (expected)
costs or invited sample size.

### Fixed total sample size

If the aim is to allocate a fixed total (invited) sample size, the basic
usage is:

    opt_nh_nonresp(
      N_h,
      phibar_h,
      S_h = NULL,
      n_total,
      ...
    )

where vectors `N_h` and `phibar_h`denote the strata population sizes,
anticipated response rates, and (optionally) strata variances,
respectively, and where scalar `n_total` denotes the total sample size.
If `S_h` is omitted, strata variances are assumed constant across
strata.

### Fixed total costs

If the aim is to allocate sample subject to a constraint on total costs,
the basic usage is:

    opt_nh_nonresp(
      N_h,
      phibar_h,
      S_h = NULL,
      cost_total,
      c_NR_h,
      tau_h,
      ...
    )

Here, `cost_total` denotes the total allowable costs, `c_NR_h` denotes
the unit costs per nonrespondent (by strata), and `tau_h` denotes the
ratio of the unit costs per respondent to those of nonrespondents (by
strata). The arguments `c_NR_h` and `tau_h` can be specified as vectors
of dimension `H` if these quantities vary by strata; alternatively, if
assumed constant across strata, they can be specified as scalars.

# References

Mendelson, J., & Elliott, M. R. (in press). Optimal allocation under
anticipated nonresponse. *Journal of Survey Statistics and Methodology*.
