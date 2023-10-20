
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

## Example

The example below shows how to use `samplingNR` to allocate a sample of
50,000 invitees, which implicitly assumes use of the constant cost per
invitee scenario (as defined in our paper) and assuming equivalent
strata variances across strata (i.e., $S_h = S$ for some constant $S$).
The [tibble](https://tibble.tidyverse.org/) `pevs_adm_2016_rrs`,
included as an example data set, contains estimated population sizes and
response propensities by strata for the 2016 Post-Election Voting Survey
of Active Duty Military (PEVS-ADM), as was used in our paper.

We can use the function `opt_nh_nonresp()` to compute the exact (i.e.,
iteratively computed) version of the optimal allocation under
nonresponse, as provided in our paper. We omit the optional argument
`S_h`, which defaults to assuming that $S_h$ is constant across strata.

``` r
library(samplingNR)
library(magrittr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

#Compute exact optimum allocation for PEVS-ADM 2016 data set for n = 50k
#Assumes tau = 1 by default
pevs_optE_alloc_50k <- opt_nh_nonresp(N_h = pevs_adm_2016_rrs$Nhat_h,
                                      phibar_h = pevs_adm_2016_rrs$rr_h,
                                      n_max = 50000)

#Merge results into original data frame
pevs_adm_alloc_50k_merged <- pevs_adm_2016_rrs %>%
  dplyr::mutate(n_h_optE = c(pevs_optE_alloc_50k),
                zeta_h_optE = attr(pevs_optE_alloc_50k,"zeta_h_min_1"))

#View results
pevs_adm_alloc_50k_merged
#> # A tibble: 91 × 10
#>    h     service paygrade age   region  sex   Nhat_h   rr_h n_h_optE zeta_h_optE
#>    <chr> <fct>   <fct>    <fct> <fct>   <fct>  <dbl>  <dbl>    <dbl>       <dbl>
#>  1 1     Army    E1-E5    18-24 US      M     1.09e5 0.0189    8267.        1.01
#>  2 2     Army    E1-E5    18-24 US      F     1.89e4 0.0364    1041.        1.03
#>  3 3     Army    E1-E5    18-24 Overse… M     1.85e4 0.0214    1333.        1.04
#>  4 4     Army    E1-E5    18-24 Overse… F     3.31e3 0.0431     179.        1.17
#>  5 5     Army    E1-E5    25-29 US      M     4.83e4 0.0487    2280.        1.01
#>  6 6     Army    E1-E5    25-29 US      F     8.68e3 0.0619     370.        1.05
#>  7 7     Army    E1-E5    25-29 Overse… M/F   1.08e4 0.0566     477.        1.04
#>  8 8     Army    E1-E5    30-34 US      M/F   2.14e4 0.0891     746.        1.01
#>  9 9     Army    E1-E5    30-34 Overse… M/F   3.74e3 0.0904     135.        1.09
#> 10 10    Army    E1-E5    35+   US      M/F   9.06e3 0.188      219.        1.02
#> # ℹ 81 more rows
```

The approximate version of this allocation (i.e., assuming that
$\zeta_h(.) \approx 1$) can be computed using the function
`opt_nh_nonresp_oneiter()`.

``` r
pevs_optA_alloc_50k <- opt_nh_nonresp_oneiter(N_h = pevs_adm_2016_rrs$Nhat_h,
                                              phibar_h = pevs_adm_2016_rrs$rr_h,
                                              n_max = 50000)
```

For comparison purposes, we can also compute Neyman allocation of
invitees ($n_h^{\textrm{Ninv}}\propto N_h S_h$), which does not adjust
for nonresponse, and Neyman allocation of expected respondents
($n_h^{\textrm{Ninv}}\propto N_h S_h / \bar{\phi}_h$), which inflates by
the inverse of the anticipated response rates. In each case, we ignore
the $S_h$ terms in our applications of these equations given our earlier
assumption that $S_h = S$ for some $S$.

The results can be seen to match the tables shown in our paper’s
Appendix.

``` r
nh_Ninv_propto <- pevs_adm_2016_rrs$Nhat_h
nh_Ninv <- nh_Ninv_propto / sum(nh_Ninv_propto) * 50000

nh_Nresp_propto <- pevs_adm_2016_rrs$Nhat_h / pevs_adm_2016_rrs$rr_h
nh_Nresp <- nh_Nresp_propto / sum(nh_Nresp_propto) * 50000

#Make sure allocations are valid
all(nh_Ninv >= 2 & nh_Ninv <= pevs_adm_2016_rrs$Nhat_h)
#> [1] TRUE
all(nh_Ninv >= 2 & nh_Nresp <= pevs_adm_2016_rrs$Nhat_h)
#> [1] TRUE

#Show results
(pevs_adm_allocs_50k <- 
  pevs_adm_alloc_50k_merged %>%
  mutate(n_h_optA = pevs_optA_alloc_50k,
         n_h_Ninv = nh_Ninv,
         n_h_Nresp = nh_Nresp) %>%
  select(c("h", "n_h_Ninv", "n_h_Nresp", "n_h_optA", "n_h_optE", "zeta_h_optE")))
#> # A tibble: 91 × 6
#>    h     n_h_Ninv n_h_Nresp n_h_optA n_h_optE zeta_h_optE
#>    <chr>    <dbl>     <dbl>    <dbl>    <dbl>       <dbl>
#>  1 1        4269.   12975.     8309.    8267.        1.01
#>  2 2         740.    1164.     1036.    1041.        1.03
#>  3 3         722.    1935.     1320.    1333.        1.04
#>  4 4         129.     172.      167.     179.        1.17
#>  5 5        1890.    2225.     2289.    2280.        1.01
#>  6 6         339.     314.      364.     370.        1.05
#>  7 7         421.     426.      473.     477.        1.04
#>  8 8         835.     537.      747.     746.        1.01
#>  9 9         146.      92.6     130.     135.        1.09
#> 10 10        354.     108.      218.     219.        1.02
#> # ℹ 81 more rows
```

Further, we can roughly examine performance by computing the number of
respondents, Kish’s design effect from weighting, and effective number
of respondents under each allocation, replicating Table 2 of our paper.
The results show that under the given assumptions, the proposed
allocation has an effective sample size that is 25% higher than the
Neyman-type allocations.

``` r
#Calculate Kish's deff_w for poststratified estimator under nonresponse
calculate_deff_w <- function (N_h, r_h) {
  stopifnot(length(N_h) == length(r_h))
  stopifnot(all(r_h >= 1))
  stopifnot(all(r_h <= N_h))
  testthat::expect_equal(all(r_h <= N_h), TRUE)

  w_h <- N_h/r_h

  wbar <- sum(w_h * r_h)/sum(r_h)
  deff_w <- 1 + 1/sum(r_h) * sum((w_h - wbar)^2 * r_h)/wbar^2

  (deff_w)
}

#Summarize results from allocation under assumption that r_h = n_h * phibar_h
summarize_alloc <- function(N_h, phibar_h, n_h) {
  n_total <- sum(n_h)
  r_total <- sum(phibar_h * n_h)
  rr <- r_total / n_total
  deff_w <- calculate_deff_w(N_h = N_h, r_h = n_h * phibar_h)
  eff_resp <- r_total / deff_w
  c(n_total, r_total, rr, deff_w, eff_resp)
}

#Summarize results of each allocation
pevs_adm_alloc_summary_list <-
  lapply(2:5,
         function(col) summarize_alloc(N_h = pevs_adm_2016_rrs$Nhat_h,
                                       phibar_h = pevs_adm_2016_rrs$rr_h,
                                       n_h = unname(unlist(pevs_adm_allocs_50k[,col]))))

names(pevs_adm_alloc_summary_list) <- names(pevs_adm_allocs_50k[,2:5])

pevs_adm_alloc_summary_table <- pevs_adm_alloc_summary_list %>%
  bind_rows() %>%
  t()

colnames(pevs_adm_alloc_summary_table) <- c("n_total", "r_total", "rr", "deff_w", "eff_resp")
tibble::as_tibble(pevs_adm_alloc_summary_table)
#> # A tibble: 4 × 5
#>   n_total r_total     rr deff_w eff_resp
#>     <dbl>   <dbl>  <dbl>  <dbl>    <dbl>
#> 1   50000   6794. 0.136    2.37    2865.
#> 2   50000   2865. 0.0573   1       2865.
#> 3   50000   4494. 0.0899   1.26    3570.
#> 4   50000   4496. 0.0899   1.26    3570.
```
