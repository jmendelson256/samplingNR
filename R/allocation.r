# Exported functions ====================

#'Compute variance inflation factor from uncertainty in responding sample size
#' (with smoothing for continuous \eqn{n_h})
#'
#'@description
#'Computes
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)},
#' as defined in our paper (Mendelson & Elliott, in press); see *Details* for a summary.
#'If there are any strata where the allocation may lead to
#' \eqn{\mathrm{E}(r_h) < r_h^{LB}} for user-specified \eqn{r_h^{LB}}
#' (3.5 by default), then \eqn{\zeta_h(.)} is evaluated
#'  using \eqn{n_h':= max\left(n_h, \left\lceil
#'      \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#' in place of \eqn{n_h}.
#'Further, \eqn{\zeta_h(.)} is computed for continuous \eqn{n_h}
#' as a weighted average
#' of evaluations at \eqn{\lfloor n_h \rfloor} and
#' \eqn{\lfloor n_h \rfloor + 1}, as in our paper.
#'
#'@details
#'In Mendelson & Elliott (in press), we assumed that the number of respondents in stratum \eqn{h}
#' can be modeled as standard binomial with support for zero removed
#' (i.e., zero-truncated binomial; see [dtruncbinom()]),
#' written as \eqn{r_h \sim TBinom(n_h, \bar{\phi}_h)}, where
#'\eqn{n_h} is the number of invitees in stratum \eqn{h},
#'\eqn{\bar{\phi}_h} is the average response propensity within stratum \eqn{h},
#'and where the unit-level response propensities are assumed constant within
#'strata.
#'Our paper defines the function
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)}.
#'This quantity is a variance inflation term that captures the effect of
#'variability in the number of respondents for a given allocation
#'(when computing the variance of the poststratified estimator under nonresponse
#'for the finite population mean).
#'
#'For discrete \eqn{n_h}, the current function (\code{calc_zeta()})
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
#'@param n_h (vector) strata sample sizes before nonresponse (\eqn{n_h})
#'@param phibar_h (vector) strata response propensities (\eqn{\bar{\phi}_h})
#'@param rh_min (scalar) minimum target respondents per stratum (\eqn{r_h^{LB}}); default is 3.5
#'@param verbose_flag (bool) flag on whether to provide noisy results
#'
#'@returns vector of length \eqn{H} containing
#'\eqn{\left\{\zeta_h(n_h', \bar{\phi}_h):h=1,2,...,H\right\}},
#'where \eqn{n_h'} is the larger of \eqn{n_h} or \eqn{\frac{r_h^{LB}}{\bar{\phi}_h}}
#'
#'@examples
#'#Basic example
#'#Note that n_h is adjusted in strata 1 and 4 since n_h * phibar_h < 3.5
#'calc_zeta(n_h = c(100, 200, 300, 300),
#'          phibar_h = c(.03, .02, .05, .005))
#'@keywords distributions
#'@references
#'Mendelson, J., & Elliott, M. R. (in press).
#  'Optimal allocation under anticipated nonresponse.
#'  \emph{Journal of Survey Statistics and Methodology}.
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

#'Calculate optimal allocation under anticipated nonresponse
#'
#'@description
#'Computes the optimal stratified sampling
#'allocation under anticipated nonresponse, as proposed in Mendelson & Elliott (in press),
#'which is
#'\deqn{n_h \propto \frac{N_h S_h \sqrt{\zeta_h(n_h, \bar{\phi}_h)}}
#'      {\sqrt{\bar{\phi}_h c_h }},}
#'where \eqn{N_h} is the stratum \eqn{h} population size,
#'\eqn{S_h} is the stratum \eqn{h} unit standard deviation,
#'\eqn{\bar{\phi}_h} is the stratum \eqn{h} average response propensity,
#'\eqn{\zeta_h(.)} is a variance inflation term that captures variability
#'in the number of respondents (see [calc_zeta()]),
#'and \eqn{c_h = \textrm{E}(C_h) / n_h = c_{NR_h} \left(\bar{\phi}_h (\tau_h - 1) + 1\right)}
#'is the expected cost per invitee in stratum \eqn{h}.
#'The cost structure assumes that respondents and nonrespondents in stratum \eqn{h}
#'have per-unit costs of \eqn{c_{R_h}} and \eqn{c_{NR_h}}, respectively,
#'with a ratio of \eqn{\tau_h = c_{R_h}/c_{NR_h}}.
#'
#'\code{opt_nh_nonresp()} computes the exact allocation in an iterative fashion
#' (see *Details* section); individual iterations are computed using
#' \code{opt_nh_nonresp_oneiter()}, which conditions on some
#' given \eqn{\zeta_h(.)}.
#'
#'Users must specify either the total sample size (\code{n_total}) or the
#' total expected costs (\code{cost_total}), but not both.
#'Under the cost specification, users must also specify the unit costs
#' of nonrespondents (\code{c_NR_h}) and the ratio of unit costs
#' for respondents to that of nonrespondents (\code{tau_h}).
#'These variables are \eqn{h}-dimensional vectors, but can be supplied as
#' scalars if equivalent across strata.
#'
#'Optionally, the user can also specify the unit population standard deviations,
#'\code{S_h}, assumed constant across strata by default.
#'
#'@details
#'
#'\code{opt_nh_nonresp()} computes the optimal allocation iteratively,
#'as follows: \enumerate{
#'  \item Iteration \eqn{k=1} calculates \eqn{n_h^1}, under the
#'  assumption that \eqn{\zeta_h(n_h, \bar{\phi}_h)=1}.
#'  \item Each subsequent iteration \eqn{k}, for \eqn{k = 2, 3, ...,}
#'  does the following: \itemize{
#'    \item Compute \eqn{\zeta_h(n_h^{k-1}, \bar{\phi}_h)} via [calc_zeta()].
#'    \item Compute \eqn{n_h^k} under the assumption that
#'    \eqn{\zeta_h(n_h, \bar{\phi}_h) =
#'    \zeta_h(n_h^{k-1}, \bar{\phi}_h)}.
#'    \item Accept the solution if the largest component of
#'     \eqn{n_h^k - n_h^{k-1}} has magnitude below some tolerance (\code{tol})
#'     or if the maximum number of iterations (\code{max_iter}) has been
#'     reached. Otherwise, continue.
#'  }
#'}
#'
#'@param N_h (vector) strata population counts (\eqn{N_h})
#'@param phibar_h (vector) strata response propensities (\eqn{\bar{\phi}_h})
#'@param ... other arguments to pass on to \code{opt_nh_nonresp_oneiter()}
#'           for applying individual iterations
#'@param tol (scalar) tolerance (for stopping)
#'@param max_iter (scalar) maximum number of iterations (>=2)
#'@param verbose_flag (bool) whether to provide detailed results
##@param browser_flag (bool) whether to open a browser for debugging purposes
#'@returns  \code{opt_nh_nonresp()} returns sample allocation
#'  vector \code{n_h} with the following attributes: \itemize{
#'  \item \code{num_iter} (scalar) number of iterations used;
#'  \item \code{zeta_h_prev} (vector) final values of \code{zeta_h} used
#'                            (i.e., from 2nd-to-last iteration); and
#'  \item \code{max_nh_delta} (scalar) biggest change in stratum
#'                             allocation from previous round.
#'}
#'
#'@examples
#' #Compute exact optimum allocation for PEVS-ADM 2016 data set for n = 50k
#' #Assumes tau = 1 by default
#' pevs_optE_alloc_50k <- opt_nh_nonresp(N_h = pevs_adm_2016_rrs$Nhat_h,
#'                                       phibar_h = pevs_adm_2016_rrs$rr_h,
#'                                       n_total = 50000)
#'
#' #Merge results into data frame
#' (pevs_adm_alloc_50k_merged <- pevs_adm_2016_rrs %>%
#'   dplyr::mutate(n_h_optE = c(pevs_optE_alloc_50k),
#'                 zeta_h_optE = attr(pevs_optE_alloc_50k,"zeta_h_prev")))
#'
#'
#'
#' #For comparison purposes, compute approximate version of proposed allocation
#' pevs_optA_alloc_50k <- opt_nh_nonresp_oneiter(N_h = pevs_adm_2016_rrs$Nhat_h,
#'                                               phibar_h = pevs_adm_2016_rrs$rr_h,
#'                                               n_total = 50000)
#'
#' #Merge results to previous tibble and reorder columns to be adjacent
#' pevs_adm_alloc_50k_merged %>%
#'   dplyr::mutate(n_h_optA = c(pevs_optA_alloc_50k)) %>%
#'   dplyr::relocate(zeta_h_optE, .after = "n_h_optA")
#'
#'@references
#'Mendelson, J., & Elliott, M. R. (in press).
#  'Optimal allocation under anticipated nonresponse.
#'  \emph{Journal of Survey Statistics and Methodology}.
#'@describeIn opt_nh_nonresp computes the exact (i.e., iterative) version of the proposed allocation.
#'@inheritDotParams opt_nh_nonresp_oneiter -zeta_h
#'@keywords allocation
#'@export
opt_nh_nonresp <- function(N_h,
                           phibar_h,
                           ...,
                           tol = 1e-8,
                           max_iter = 20,
                           verbose_flag = FALSE) {
  #if(browser_flag) browser()
  stopifnot(length(max_iter)==1)
  stopifnot(max_iter>=2)

  dotdotdot_var_in_formals <- names(list(...)) %in% names(formals(opt_nh_nonresp_oneiter))
  if(!all(dotdotdot_var_in_formals)) {
    vars_not_in_formals <- paste0(names(list(...))[!dotdotdot_var_in_formals], collaose = ", ")
    stop(paste0(c("Contains extraneous vars not in opt_nh_nonresp_oneiter: ", vars_not_in_formals),
                collapse = ""))
  }
  if("zeta_h" %in% names(list(...))) stop("zeta should be NULL; will be calculated in procedure")

  nh_list <- vector(mode = "list", length = max_iter)
  zeta_list <- vector(mode = "list", length = max_iter)

  #Compute initial allocation
  nh_args <- append(list(N_h = N_h,
                         phibar_h = phibar_h),
                    list(...))

  nh_list[[1]] <-  do.call(opt_nh_nonresp_oneiter, nh_args)
  zeta_list[[1]] <- calc_zeta(n_h = nh_list[[1]],
                              phibar_h = phibar_h)
  iter <- 1
  max_nh_delta <- Inf

  while(iter < max_iter & max_nh_delta > tol) {
    iter <- iter + 1

    nh_args <- append(list(N_h = N_h,
                           phibar_h = phibar_h,
                           zeta_h = zeta_list[[iter - 1]]),
                      list(...))
    nh_list[[iter]] <-  do.call(opt_nh_nonresp_oneiter, nh_args)

    zeta_list[[iter]] <- calc_zeta(n_h = nh_list[[iter]],
                                   phibar_h = phibar_h)

    max_nh_delta <- max(abs(nh_list[[iter]] - nh_list[[iter-1]]))
  }

  nh_list <- nh_list[1:iter]
  zeta_list <- zeta_list[1:iter]

  n_h <- nh_list[[iter]]
  zeta_h_prev <- zeta_list[[iter-1]]

  res <- n_h
  attr(res, "num_iter") <- iter
  attr(res, "zeta_h_prev") <- zeta_h_prev
  attr(res, "max_nh_delta") <- max_nh_delta

  if(verbose_flag) {
    return(mget(ls()))
  } else {
    return(res)
  }
}

# Internal functions ====================

#'Compute variance inflation factor from uncertainty in responding sample size
#' (for integer \eqn{n_h})
#'
#'@description
#'Computes
#'\eqn{\zeta_h(n_h, \bar{\phi}_h):=\mathrm{E}(r_h)\mathrm{E}\left(\frac{1}{r_h}\right)},
#' as defined in the paper.
#'If there are any strata where the allocation may lead to
#' \eqn{\mathrm{E}(r_h) < r_h^{LB}} for user-specified \eqn{r_h^{LB}}
#' (3.5 by default), then \eqn{\zeta_h(.)} is evaluated
#'  using \eqn{n_h':= max\left(n_h, \left\lceil
#'      \frac{r_h^{LB}}{\bar{\phi}_h}\right\rceil\right)}
#' in place of \eqn{n_h}.
#'Continuous values of \eqn{n_h} are rounded (by default).
#'
#'@details
#'In our paper, we assumed that the number of respondents in stratum \eqn{h}
#' can be modeled as standard binomial with support for zero removed
#' (i.e., zero-truncated binomial; see [dtruncbinom()]),
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
#'The current function (\code{calc_zeta_discrete}())
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
#'@param n_h (vector) strata sample sizes before nonresponse (\eqn{n_h})
#'@param phibar_h (vector) strata response propensities (\eqn{\bar{\phi}_h})
#'@param rh_min (scalar) minimum target respondents per stratum (\eqn{r_h^{LB}}); default is 3.5
#'@param round_flag (bool) specifies whether to round continuous allocations
#'                         (versus throwing an error)
#'@param verbose_flag (bool) flag on whether to provide detailed results
#'
#'@returns vector of length \eqn{H} containing
#'\eqn{\left\{\zeta_h(n_h', \bar{\phi}_h):h=1,2,...,H\right\}},
#'where \eqn{n_h'} is the larger of \eqn{n_h} or \eqn{\frac{r_h^{LB}}{\bar{\phi}_h}}
#'
#'@examples
#'\dontrun{
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
#'calc_zeta_discrete(n_h = c(100, 200, 300, 300),
#'                   phibar_h = c(.03, .02, .05, .001),
#'                   rh_min = 0,
#'                   verbose_flag = TRUE)
#'}
#'@seealso [calc_zeta()] for use with continuous \code{n_h}
#'
#'@keywords internal distributions
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

#'Calculates one iteration of proposed optimal allocation conditioned on some zeta
#'
#'@description
#'\code{opt_nh_nonresp_oneiter()} computes one iteration of the proposed
#'allocation for some user-supplied \eqn{\zeta_h(.)}.
#'
#'@param N_h (vector) strata population counts (\eqn{N_h})
#'@param phibar_h (vector) strata response propensities (\eqn{\bar{\phi}_h})
#'@param S_h (vector) strata population standard deviations (\eqn{S_h}); constant, by default
#'@param n_total (scalar) total sample to allocate
#'@param cost_total (scalar) total expected costs to allocate
#'@param zeta_h (vector; use with \code{opt_nh_nonresp_oneiter}, only;
#'           optional) adjustment factor to reflect inflation in variances
#'           from randomness in the number of respondents (default = 1)
#'@param c_NR_h (vector; use with cost_total) per-unit costs for nonrespondents in stratum h
#'  (\eqn{c_{NR_h}})
#'@param tau_h (vector; use with cost_total) ratio of costs for respondents to
#'  costs for nonrespondents in stratum h (\eqn{\tau_h})
#'@param strict_flag (boolean) whether to throw error (versus warning) if any
#'  \eqn{n_h > N_h}
#'@param verbose_flag (boolean) whether to provide detailed results
#'
#Old examples (not currently used):
#
# #Basic example with varying response rates
# opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
# phibar = c(.2, .1, .05),
# n_total = 5000)
#
# #Expanded example with more varying factors
# opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
#                           phibar_h = c(.2, .1, .05),
#                           tau_h = c(1.1, 1.2, 1.3),
#                           c_NR_h = 2:4,
#                           n_total = 5e4)
#
# #Same as above but fix total costs
# opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
#                           phibar_h = c(.2, .1, .05),
#                           tau_h = c(1.1, 1.2, 1.3),
#                           c_NR_h = 2:4,
#                           cost_total = 5e4)
#
# \dontrun{
# #gives warning due to excessive n_h's
# opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
#                           phibar_h = c(.2, .1, .05),
#                           tau_h = c(1.1, 1.2, 1.3),
#                           c_NR_h = 2:4,
#                           cost_total = 5e6)
# }
#'@returns \code{opt_nh_nonresp_oneiter} returns sample allocation vector
#'         \code{n_h}, computed from a single iteration given
#'         the user-supplied \code{zeta_h}.
#'
#'@describeIn opt_nh_nonresp computes a single iteration of the
#'            proposed allocation for user-supplied \code{zeta_h}
#'            (used in place of the \eqn{\zeta_h(.)} term).
#'@keywords allocation
#'@export
opt_nh_nonresp_oneiter <- function(N_h,
                                   phibar_h,
                                   S_h = NULL,
                                   n_total = NULL,
                                   cost_total = NULL,
                                   zeta_h = NULL,
                                   c_NR_h = NULL,
                                   tau_h = NULL,
                                   strict_flag = TRUE,
                                   verbose_flag = FALSE) {

  #Check object lengths
  H <- length(N_h)
  if(length(N_h) != length(phibar_h)) stop("N_h and phibar_h are different lengths")

  if(!xor(!is.null(n_total),
          !is.null(cost_total))) {
    stop("User must specify either n_total or cost_total (and not both)")
  }

  if(is.null(n_total)) {
    #If allocating based on costs, user must input costs
    if(is.null(c_NR_h)) stop("User must specify n_total or provide cost information (c_NR_h is missing)")
    if(is.null(tau_h)) stop("User must specify n_total or provide cost information (tau_h is missing)")

    #if c_NR_h or tau_h are scalar, apply to all strata
    if(length(c_NR_h)==1) c_NR_h <- rep(c_NR_h, H)
    if(length(tau_h)==1) tau_h <- rep(tau_h, H)
  } else {
    if(!is.null(tau_h)) stop("tau_h is only intended for use with cost_total, but n_total was provided")
    if(!is.null(c_NR_h)) stop("c_NR_h is only intended for use with cost_total, but n_total was provided")
  }

  if(is.null(S_h)) S_h <- rep(1, H)
  if(is.null(zeta_h)) zeta_h <- rep(1, H)
  if(is.null(c_NR_h)) c_NR_h <- rep(1, H)
  if(is.null(tau_h)) tau_h <- rep(1, H)
  if(length(N_h) != length(S_h)) {
    stop("N_h and S_h are different lengths")
  }
  if(length(N_h) != length(zeta_h)) {
    stop("N_h and zeta_h are different lengths")
  }
  if(length(N_h) != length(c_NR_h)) {
    stop("N_h and c_NR_h are different lengths")
  }
  if(length(N_h) != length(tau_h)) {
    stop("N_h and tau_h are different lengths")
  }

  #Check for illogical arguments
  if(any(phibar_h <= 0)) stop("invalid phibar_h (must be strictly positive)")
  if(any(phibar_h > 1)) stop("invalid phibar_h (cannot exceed 1)")
  if(any(S_h <= 0)) stop("invalid S_h (must be strictly positive)")
  if(any(c_NR_h <= 0)) stop("invalid c_NR_h (must be strictly positive)")
  if(any(tau_h <= 0)) stop("invalid tau_h (must be strictly positive)")

  nh_propto_num <- N_h * S_h * sqrt(zeta_h)
  nh_propto_denom <- sqrt(phibar_h * c_NR_h * (phibar_h * tau_h + 1 - phibar_h))
  nh_propto <- nh_propto_num / nh_propto_denom

  nh_one_invitee <- nh_propto / sum(nh_propto) #share of sample for strat h

  if(!is.null(n_total)) {
    n_h <- nh_one_invitee * n_total
    if(!exists("C_h")) C_h <- NULL #allows verbose results even if costs not explicitly specified
  } else {
    exp_cost_per_invitee <- sum(nh_one_invitee * c_NR_h * (phibar_h * tau_h + 1 - phibar_h))
    n_h <- nh_one_invitee * cost_total / exp_cost_per_invitee
    C_h <- n_h * c_NR_h * (phibar_h * tau_h + 1 - phibar_h)
  }

  res <- n_h

  if(!all(n_h <= N_h)) {
    my_error_msg <- "some n_h's exceed N_h's; try mathematical programming approach"
    if(strict_flag) {
      stop(my_error_msg)
    } else {
      warning(my_error_msg)
    }
  }

  if(verbose_flag) {
    results_tbl <- tibble::tibble(N_h = N_h,
                                  S_h = S_h,
                                  zeta_h = zeta_h,
                                  c_NR_h = c_NR_h,
                                  tau_h = tau_h,
                                  phibar_h = phibar_h,
                                  n_h = n_h,
                                  C_h = C_h)
    return(mget(ls()))
  } else {
    return(n_h)
  }
}
#codetools::findGlobals(opt_nh_nonresp_oneiter)
