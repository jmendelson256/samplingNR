% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncbinom.R
\name{dtruncbinom}
\alias{dtruncbinom}
\title{Calculate density of truncated binomial distribution}
\usage{
dtruncbinom(x, size, prob)
}
\arguments{
\item{x}{vector of values}

\item{size}{(scalar) number of trials (\code{n}); must be 1 or more}

\item{prob}{(scalar) probability of success on each trial (\code{p}); must be nonzero}
}
\value{
Vector of densities associated with the provided values (x) for given parameters
}
\description{
Calculates density of the binomial distribution where support for
0 has been removed, and which we refer to as the
\emph{(zero-) truncated binomial distribution} (e.g., Rider, 1955;
Stephan, 1945).
}
\details{
Let \eqn{X} denote a random variable with probability mass function
(pmf) of

\eqn{\mathrm{p}(k;n,p) = \mathrm{Pr}(X = k; n, p) =
\frac{\binom{n}{k} p^k (1-p)^{n-k}}{1 - (1-p)^n} \propto
\binom{n}{k} p^k (1-p)^{n-k}},

for \eqn{k = 1, 2, ..., n}, and with zero mass, otherwise.
We say that \eqn{X} is \emph{truncated binomial} with parameters \eqn{(n,p)},
written as \eqn{X \sim TBinom(n,p)},
and where \eqn{n} refers to the number of trials and
\eqn{p} refers to the probability of success of each trial for the
corresponding binomial distribution.
}
\examples{
dtruncbinom(1:20, 20, .2)

#same as above
dbinom(1:20, size = 20, prob = .2) / sum(dbinom(1:20, size = 20, prob = .2))

}
\references{
Rider, P. R. (1955). Truncated binomial and negative binomial distributions.
\emph{Journal of the American Statistical Association, 50}(271), 877-883.

Stephan, F. F. (1945). The expected value and variance of the reciprocal and
other negative powers of a positive Bernoullian variate.
\emph{The Annals of Mathematical Statistics, 16}(1), 50-61.
}
\seealso{
\code{\link[=calc_moments_truncbinom]{calc_moments_truncbinom()}} for computing TBinom moments
}
\keyword{distributions}
