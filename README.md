
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplingNR

<!-- badges: start -->
<!-- badges: end -->

`samplingNR` is an R package that allows for computing optimal
allocations under anticipated nonresponse.

## Installation

You can install the development version of `samplingNR` by using a
command such as:

``` r
#install.packages("samplingNR_0.2.0.tar.gz")
```

## Vignette

You can learn about how to use the package in `vignette("samplingNR")`,
which shows how to replicate two of the tables from our paper in an
application to a post-election survey of military personnel.

## Basic usage

The main allocation function is `opt_nh_nonresp()`, which provides the
exact version of our proposed optimal allocation under anticipated
nonresponse. This function provides the optimal allocation subject to
constraints on total sample size or on total costs.

### Fixed total sample size

If the aim is to allocate a fixed total (invited) sample size, the basic
usage is:

    opt_nh_nonresp(
      N_h,
      phibar_h,
      S_h = NULL,
      n_max
      ...
    )

where vectors `N_h`, `phibar_h`, and `S_h` denote the strata population
sizes, anticipated response rates, and (optionally) strata variances,
respectively, and where `n_max` denotes the total sample size. If `S_h`
is omitted, strata variances are assumed constant across strata.

### Fixed total costs

If the aim is to allocate sample subject to a constraint on total costs,
the basic usage is:

    opt_nh_nonresp(
      N_h,
      phibar_h,
      S_h = NULL,
      c_max,
      c_NR_h,
      tau_h
      ...
    )

Here, `c_max` denotes the total allowable costs, `c_NR_h` denotes the
unit costs per nonrespondent (by strata), and `tau_h` denotes the ratio
of the unit costs per respondent to those of nonrespondents (by strata).
The arguments `c_NR_h` and `tau_h` can be specified as vectors of
dimension `H` if these quantities vary by strata; alternatively, if
assumed constant across strata, they can be specified as scalar.
