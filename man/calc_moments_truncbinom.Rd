% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncbinom.R
\name{calc_moments_truncbinom}
\alias{calc_moments_truncbinom}
\title{Compute moments for truncated binomial distribution}
\usage{
calc_moments_truncbinom(size, prob)
}
\arguments{
\item{size}{scalar number of trials (\code{n}); must be 1 or more}

\item{prob}{scalar probability of success on each trial (\code{p}); must be nonzero}
}
\value{
vector containing \eqn{\textrm{E}(X)}, \eqn{\textrm{E}(X^2)},
\eqn{\textrm{Var}(X)}, and \eqn{\textrm{E}(1/x)}
}
\description{
For random variable \eqn{X \sim TBinom(n,p)}, which follows the
(zero)-truncated binomial distribution as defined
by \code{\link[=dtruncbinom]{dtruncbinom()}} (i.e., binomial
distribution with support for zero removed),
\code{\link[=calc_moments_truncbinom]{calc_moments_truncbinom()}} computes \eqn{\textrm{E}(X)}, \eqn{\textrm{E}(X^2)},
\eqn{\textrm{Var}(X)}, and \eqn{\textrm{E}(1/X)}.
}
\details{
Note that moments are computed directly using the probability mass function since there
is not a simple closed form for \eqn{\textrm{E}(1/x)}.
This may lead to slow computation for extremely large numbers of trials.
}
\examples{
calc_moments_truncbinom(30, .1)
calc_moments_truncbinom(30, .2)

}
\seealso{
\code{\link[=dtruncbinom]{dtruncbinom()}} for computing TBinom densities
}
\keyword{distributions}
