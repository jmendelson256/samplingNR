#'Compute zeta, as defined in paper, with smoothing for fractional n
#'
#'@description
#'\code{calc_zeta} computes
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)},
#' as defined in the paper.
#'If there are any strata where the allocation may lead to
#' \eqn{\mathrm{E}(r_h) < 3.5}, we evaluate
#' \eqn{\zeta_h(.)} using
#' \eqn{n_h':= max\left(n_h, \left\lceil
#'      \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#' in place of \eqn{n_h}.
#'Further, \eqn{\zeta_h(.)} is computed for continuous \eqn{n_h}
#' as a weighted average
#' of evaluations at \eqn{\lfloor n_h \rfloor} and
#' \eqn{\lfloor n_h \rfloor + 1}, as in the paper.
#'
#'@details
#'In our paper, we assumed that the number of respondents in stratum \eqn{h}
#' can be modeled as standard binomial with support for zero removed
#' (i.e., zero-truncated binomial; see \code{\link{dtruncbinom}}),
#' written as \eqn{r_h \sim TBinom(n_h, \bar{\phi}_h)}, where
#'\eqn{n_h} is the number of invitees in stratum \eqn{h},
#'\eqn{\bar{\phi}_h} is the average response propensity within stratum \eqn{h},
#'and where the unit-level response propensities are assumed constant within
#'strata.
#'Our paper defines the function
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)}.
#'This quantity is a variance inflation term that captures the effect of
#'variability in the number of respondents for a given allocation.
#'
#'For discrete \eqn{n_h}, the current function (\code{calc_zeta})
#'computes \eqn{\zeta_h(n_h', \bar{\phi}_h)}, where we use
#'\eqn{n_h':= max\left(n_h, \left\lceil \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#'to avoid underallocating to strata with too few expected respondents, and
#'where \eqn{r_h^{LB}} is some given lower bound on the number
#'of expected respondents.
#'By default, we set \eqn{r_h^{LB} = 3.5}, since the truncated binomial
#'distribution may sometimes be a poor approximation
#'for the binomial distribution below these levels,
#'and as we observed numerically that \eqn{\zeta_h(n_h, \bar{\phi}_h)}
#'is roughly maximized for various \eqn{\bar{\phi}_h} (fixed at levels
#'between .01 and 1) when \eqn{n_h \approx \frac{3.5}{\bar{\phi}_h}}.
#'
#'For continuous \eqn{n_h}, we define \eqn{\zeta_h(n_h, \bar{\phi}_h)}
#' as a weighted average of its evaluations at \eqn{\lfloor n_h \rfloor} and
#' \eqn{\lfloor n_h \rfloor + 1}, via
#'\deqn{\zeta_h'(n_h,\bar{\phi}_h) =
#'     w_h \cdot \zeta_h(\lfloor n_h \rfloor,\bar{\phi}_h) + \left(1 - w_h\right)\cdot
#'     \zeta_h(\lfloor n_h \rfloor + 1,\bar{\phi}_h),}
#' where \eqn{w_h= \left(\lfloor n_h \rfloor + 1\right) - n_h}.
#'
#'@param n_h (vector) strata sample sizes (before nonresponse)
#'@param phibar_h (vector) strata response propensities
#'@param rh_min (scalar) minimum target respondents per stratum (default 3.5)
#'@param verbose_flag (bool) flag on whether to provide noisy results
#'
#'@returns vector of length \eqn{H} containing
#'\eqn{\left\{\zeta_h(n_h', \bar{\phi}_h):h=1,2,...,H\right\}},
#'where \eqn{n_h'} is the larger of \eqn{n_h} or \eqn{\frac{r_h^{LB}}{\bar{\phi}_h}}.
#'
#'@examples
#'#Basic example
#'#Note that n_h is adjusted in strata 1 and 4 since n_h * phibar_h < 3.5
#'calc_zeta(n_h = c(100, 200, 300, 300),
#'          phibar_h = c(.03, .02, .05, .005))
#'
#'@export
calc_zeta <- function(n_h,
                      phibar_h,
                      rh_min = 3.5,
                      verbose_flag = FALSE) {
  #note: n_h_ceil is n_h_floor + 1 instead of ceiling(n_h_floor) to allow for integer inputs
  n_h_floor <- floor(n_h)
  n_h_ceil <- n_h_floor + 1
  n_h_wt <- 1 - (n_h - n_h_floor)

  zeta_floor <- calc_zeta_discrete(n_h = n_h_floor,
                                   phibar_h = phibar_h,
                                   rh_min = rh_min,
                                   verbose_flag = verbose_flag)

  zeta_ceil <- calc_zeta_discrete(n_h = n_h_ceil,
                                  phibar_h = phibar_h,
                                  rh_min = rh_min,
                                  verbose_flag = verbose_flag)

  if(verbose_flag) {
    res <- n_h_wt  * zeta_floor$res + (1 - n_h_wt) * zeta_ceil$res
  } else {
    res <- n_h_wt  * zeta_floor + (1 - n_h_wt) * zeta_ceil
  }

  if(verbose_flag) {
    return(mget(ls()))
  } else {
    return(res)
  }
}
#codetools::findGlobals(calc_zeta)

#'Compute zeta, as defined in paper, for integer number of respondents
#'
#'@description
#'\code{calc_zeta_discrete} computes
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)},
#'as defined in the paper, for integer \eqn{n_h}.
#'If there are any strata where the allocation may lead to
#' \eqn{\mathrm{E}(r_h) < 3.5}, we evaluate
#' \eqn{\zeta_h(.)} using
#' \eqn{n_h':= max\left(n_h, \left\lceil
#'      \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#' in place of \eqn{n_h}.
#'Continuous values of \eqn{n_h} are rounded (by default).
#'
#'@details
#'In our paper, we assumed that the number of respondents in stratum \eqn{h}
#' can be modeled as standard binomial with support for zero removed
#' (i.e., zero-truncated binomial; see \code{\link{dtruncbinom}}),
#' written as \eqn{r_h \sim TBinom(n_h, \bar{\phi}_h)}, where
#'\eqn{n_h} is the number of invitees in stratum \eqn{h},
#'\eqn{\bar{\phi}_h} is the average response propensity within stratum \eqn{h},
#'and where the unit-level response propensities are assumed constant within
#'strata.
#'Our paper defines the function
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)}.
#'This quantity is a variance inflation term that captures the effect of
#'variability in the number of respondents for a given allocation.
#'
#'The current function (\code{calc_zeta_discrete})
#'computes \eqn{\zeta_h(n_h', \bar{\phi}_h)}, where we use
#'\eqn{n_h':= max\left(n_h, \left\lceil \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#'to avoid underallocating to strata with too few expected respondents, and
#'where \eqn{r_h^{LB}} is some given lower bound on the number
#'of expected respondents.
#'By default, we set \eqn{r_h^{LB} = 3.5}, since the truncated binomial
#'distribution may sometimes be a poor approximation
#'for the binomial distribution below these levels,
#'and as we observed numerically that \eqn{\zeta_h(n_h, \bar{\phi}_h)}
#'is roughly maximized for various \eqn{\bar{\phi}_h} (fixed at levels
#'between .01 and 1) when \eqn{n_h \approx \frac{3.5}{\bar{\phi}_h}}.
#'
#'Function is only defined for discrete \eqn{n_h}; input \code{round_flag}
#'controls whether \eqn{n_h} is rounded to the nearest integer
#'(as opposed to throwing an error).
#'
#'@param n_h (vector) strata sample sizes (before nonresponse)
#'@param phibar_h (vector) strata response propensities
#'@param rh_min (scalar) minimum target respondents per stratum (default 3.5)
#'@param round_flag (bool) specifies whether to round continuous allocations
#'                         (versus throwing an error)
#'@param verbose_flag (bool) flag on whether to provide detailed results
#'
#'@returns vector of length \eqn{H} containing
#'\eqn{\left\{\zeta_h(n_h', \bar{\phi}_h):h=1,2,...,H\right\}},
#'where \eqn{n_h'} is the larger of \eqn{n_h} or \eqn{\frac{r_h^{LB}}{\bar{\phi}_h}}.
#'
#'@examples
#'
#'#Basic example
#'calc_zeta_discrete(n_h = c(100, 200, 300, 300),
#'                   phibar_h = c(.03, .02, .05, .001))
#'
#'#Verbose info shows that strata 1,5 have their sample sizes increased
#'#   due to having fewer than 3.5 expected respondents
#' calc_zeta_discrete(n_h = c(100, 200, 300, 300),
#'                    phibar_h = c(.03, .02, .05, .001),
#'                    verbose_flag = TRUE)
#'
#'#Similar to above, but removes the minimum number of respondents,
#'# which changes the values for strata 1 and 4.
#' calc_zeta_discrete(n_h = c(100, 200, 300, 300),
#'                    phibar_h = c(.03, .02, .05, .001),
#'                    rh_min = 0,
#'                    verbose_flag = TRUE)
#'@keywords internal
calc_zeta_discrete <- function(n_h,
                               phibar_h,
                               rh_min = 3.5,
                               round_flag = TRUE,
                               verbose_flag = FALSE) {
  stopifnot(length(rh_min)==1)
  stopifnot(length(n_h)==length(phibar_h))
  stopifnot(all(phibar_h>0))
  stopifnot(all(phibar_h<=1))

  nh_orig <- n_h
  nh_min <- ceiling(rh_min / phibar_h)
  if(round_flag) {
    n_h <- round(n_h)
  } else {
    stopifnot(all(n_h == as.integer(n_h)))
  }
  n_h[nh_orig < nh_min] <- nh_min[nh_orig < nh_min]


  rh_moments <-
    lapply(1:length(n_h),
           function(i) calc_moments_truncbinom(size = n_h[i], prob = phibar_h[i])) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(h = 1:length(.env$n_h),
                  n_h = .env$n_h,
                  phibar_h = .env$phibar_h,
                  zeta_h = .data$E_X * .data$E_one_over_X) %>%
    dplyr::relocate(.data$h)

  res <- rh_moments$zeta_h

  if(verbose_flag) {
    return(mget(ls()))
  } else {
    return(res)
  }
}
#codetools::findGlobals(calc_zeta_discrete)
