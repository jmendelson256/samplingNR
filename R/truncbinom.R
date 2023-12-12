# Truncated binomial distribution functions ========

#'Calculate density of truncated binomial distribution
#'
#'@description
#'
#'Calculates density of the binomial distribution where support for
#' 0 has been removed, and which we refer to as the
#' \emph{(zero-) truncated binomial distribution} (e.g., Rider, 1955;
#' Stephan, 1945).
#'
#'@details
#'Let \eqn{X} denote a random variable with probability mass function
#' (pmf) of
#'
#'\eqn{\mathrm{p}(k;n,p) = \mathrm{Pr}(X = k; n, p) =
#' \frac{\binom{n}{k} p^k (1-p)^{n-k}}{1 - (1-p)^n} \propto
#' \binom{n}{k} p^k (1-p)^{n-k}},
#'
#'for \eqn{k = 1, 2, ..., n}, and with zero mass, otherwise.
#' We say that \eqn{X} is \emph{truncated binomial} with parameters \eqn{(n,p)},
#'  written as \eqn{X \sim TBinom(n,p)},
#'  and where \eqn{n} refers to the number of trials and
#'  \eqn{p} refers to the probability of success of each trial for the
#'  corresponding binomial distribution.
#'
#'@param x vector of values
#'@param size (scalar) number of trials (\code{n}); must be 1 or more
#'@param prob (scalar) probability of success on each trial (\code{p}); must be nonzero
#'
#'@returns Vector of densities associated with the provided values (x) for given parameters
#'
#'@examples
#'dtruncbinom(1:20, 20, .2)
#'
#'#same as above
#'dbinom(1:20, size = 20, prob = .2) / sum(dbinom(1:20, size = 20, prob = .2))
#'
#'@references
#'Rider, P. R. (1955). Truncated binomial and negative binomial distributions.
#' \emph{Journal of the American Statistical Association, 50}(271), 877-883.
#'
#'Stephan, F. F. (1945). The expected value and variance of the reciprocal and
#' other negative powers of a positive Bernoullian variate.
#' \emph{The Annals of Mathematical Statistics, 16}(1), 50-61.
#'
#'@seealso \code{\link[=calc_moments_truncbinom]{calc_moments_truncbinom()}} for computing TBinom moments
#'@keywords distributions
#'@export
dtruncbinom <- function(x, size, prob) {
  stopifnot(length(size)==1)
  stopifnot(length(prob)==1)
  stopifnot(size >= 1)
  stopifnot(prob > 0)

  res <- (stats::dbinom(x = x, size = size, prob = prob) / (1 - (1 - prob)^size))
  res[x==0] <- 0
  (res)
}


#'Compute moments for truncated binomial distribution
#'
#'For random variable \eqn{X \sim TBinom(n,p)}, which follows the
#'(zero)-truncated binomial distribution as defined
#' by [dtruncbinom()] (i.e., binomial
#' distribution with support for zero removed),
#' [calc_moments_truncbinom()] computes \eqn{\textrm{E}(X)}, \eqn{\textrm{E}(X^2)},
#' \eqn{\textrm{Var}(X)}, and \eqn{\textrm{E}(1/X)}.
#'
#'@details
#'Note that moments are computed directly using the probability mass function since there
#'is not a simple closed form for \eqn{\textrm{E}(1/x)}.
#'This may lead to slow computation for extremely large numbers of trials.
#'
#'@param size scalar number of trials (\code{n}); must be 1 or more
#'@param prob scalar probability of success on each trial (\code{p}); must be nonzero
#'
#'@returns vector containing \eqn{\textrm{E}(X)}, \eqn{\textrm{E}(X^2)},
#' \eqn{\textrm{Var}(X)}, and \eqn{\textrm{E}(1/x)}
#'
#'@examples
#'calc_moments_truncbinom(30, .1)
#'calc_moments_truncbinom(30, .2)
#'
#'@seealso \code{\link[=dtruncbinom]{dtruncbinom()}} for computing TBinom densities
#'@keywords distributions
#'@export
calc_moments_truncbinom <- function(size, prob) {
  stopifnot(length(size)==1)
  stopifnot(length(prob)==1)

  my_pdf <- function(x) dtruncbinom(x, size = size, prob = prob)

  my_density_vec <- my_pdf(1:size)

  E_X <- sum(1:size * my_density_vec)
  E_X_sq <- sum((1:size)^2 * my_density_vec)
  V_X <- sum((1:size - E_X)^2 * my_density_vec)
  E_one_over_X <- sum(my_density_vec / 1:size)

  #my_pdf <- dtruncbinom_factory(size = size, prob = prob)

  res <- c(E_X = E_X,
           E_X_sq = E_X_sq,
           V_X = V_X,
           E_one_over_X = E_one_over_X)

  (res)
}

