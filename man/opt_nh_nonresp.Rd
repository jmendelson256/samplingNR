% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allocation.r
\name{opt_nh_nonresp}
\alias{opt_nh_nonresp}
\alias{opt_nh_nonresp_oneiter}
\title{Calculate optimal allocation under anticipated nonresponse}
\usage{
opt_nh_nonresp(
  N_h,
  phibar_h,
  ...,
  tol = 1e-08,
  max_iter = 20,
  verbose_flag = FALSE
)

opt_nh_nonresp_oneiter(
  N_h,
  phibar_h,
  S_h = NULL,
  n_total = NULL,
  cost_total = NULL,
  zeta_h = NULL,
  c_NR_h = NULL,
  tau_h = NULL,
  strict_flag = TRUE,
  verbose_flag = FALSE
)
}
\arguments{
\item{N_h}{(vector) strata population counts (\eqn{N_h})}

\item{phibar_h}{(vector) strata response propensities (\eqn{\bar{\phi}_h})}

\item{...}{
  Arguments passed on to \code{\link[=opt_nh_nonresp_oneiter]{opt_nh_nonresp_oneiter}}
  \describe{
    \item{\code{}}{}
  }}

\item{tol}{(scalar) tolerance (for stopping)}

\item{max_iter}{(scalar) maximum number of iterations (>=2)}

\item{verbose_flag}{(boolean) whether to provide detailed results}

\item{S_h}{(vector) strata population standard deviations (\eqn{S_h}); constant, by default}

\item{n_total}{(scalar) total sample to allocate}

\item{cost_total}{(scalar) total expected costs to allocate}

\item{zeta_h}{(vector; use with \code{opt_nh_nonresp_oneiter}, only;
optional) adjustment factor to reflect inflation in variances
from randomness in the number of respondents (default = 1)}

\item{c_NR_h}{(vector; use with cost_total) per-unit costs for nonrespondents in stratum h
(\eqn{c_{NR_h}})}

\item{tau_h}{(vector; use with cost_total) ratio of costs for respondents to
costs for nonrespondents in stratum h (\eqn{\tau_h})}

\item{strict_flag}{(boolean) whether to throw error (versus warning) if any
\eqn{n_h > N_h}}
}
\value{
\code{opt_nh_nonresp()} returns sample allocation
vector \code{n_h} with the following attributes: \itemize{
\item \code{num_iter} (scalar) number of iterations used;
\item \code{zeta_h_prev} (vector) final values of \code{zeta_h} used
(i.e., from 2nd-to-last iteration); and
\item \code{max_nh_delta} (scalar) biggest change in stratum
allocation from previous round.
}

\code{opt_nh_nonresp_oneiter} returns sample allocation vector
\code{n_h}, computed from a single iteration given
the user-supplied \code{zeta_h}.
}
\description{
Computes the optimal stratified sampling
allocation under anticipated nonresponse, as proposed in Mendelson & Elliott (in press),
which is
\deqn{n_h \propto \frac{N_h S_h \sqrt{\zeta_h(n_h, \bar{\phi}_h)}}
     {\sqrt{\bar{\phi}_h c_h }},}
where \eqn{N_h} is the stratum \eqn{h} population size,
\eqn{S_h} is the stratum \eqn{h} unit standard deviation,
\eqn{\bar{\phi}_h} is the stratum \eqn{h} average response propensity,
\eqn{\zeta_h(.)} is a variance inflation term that captures variability
in the number of respondents (see \code{\link[=calc_zeta]{calc_zeta()}}),
and \eqn{c_h = \textrm{E}(C_h) / n_h = c_{NR_h} \left(\bar{\phi}_h (\tau_h - 1) + 1\right)}
is the expected cost per invitee in stratum \eqn{h}.
The cost structure assumes that respondents and nonrespondents in stratum \eqn{h}
have per-unit costs of \eqn{c_{R_h}} and \eqn{c_{NR_h}}, respectively,
with a ratio of \eqn{\tau_h = c_{R_h}/c_{NR_h}}.

\code{opt_nh_nonresp()} computes the exact allocation in an iterative fashion
(see \emph{Details} section); individual iterations are computed using
\code{opt_nh_nonresp_oneiter()}, which conditions on some
given \eqn{\zeta_h(.)}.

Users must specify either the total sample size (\code{n_total}) or the
total expected costs (\code{cost_total}), but not both.
Under the cost specification, users must also specify the unit costs
of nonrespondents (\code{c_NR_h}) and the ratio of unit costs
for respondents to that of nonrespondents (\code{tau_h}).
These variables are \eqn{h}-dimensional vectors, but can be supplied as
scalars if equivalent across strata.

Optionally, the user can also specify the unit population standard deviations,
\code{S_h}, assumed constant across strata by default.

\code{opt_nh_nonresp_oneiter()} computes one iteration of the proposed
allocation for some user-supplied \eqn{\zeta_h(.)}.
}
\details{
\code{opt_nh_nonresp()} computes the optimal allocation iteratively,
as follows: \enumerate{
\item Iteration \eqn{k=1} calculates \eqn{n_h^1}, under the
assumption that \eqn{\zeta_h(n_h, \bar{\phi}_h)=1}.
\item Each subsequent iteration \eqn{k}, for \eqn{k = 2, 3, ...,}
does the following: \itemize{
\item Compute \eqn{\zeta_h(n_h^{k-1}, \bar{\phi}_h)} via \code{\link[=calc_zeta]{calc_zeta()}}.
\item Compute \eqn{n_h^k} under the assumption that
\eqn{\zeta_h(n_h, \bar{\phi}_h) =
   \zeta_h(n_h^{k-1}, \bar{\phi}_h)}.
\item Accept the solution if the largest component of
\eqn{n_h^k - n_h^{k-1}} has magnitude below some tolerance (\code{tol})
or if the maximum number of iterations (\code{max_iter}) has been
reached. Otherwise, continue.
}
}
}
\section{Functions}{
\itemize{
\item \code{opt_nh_nonresp()}: computes the exact (i.e., iterative) version of the proposed allocation.

\item \code{opt_nh_nonresp_oneiter()}: computes a single iteration of the
proposed allocation for user-supplied \code{zeta_h}
(used in place of the \eqn{\zeta_h(.)} term).

}}
\examples{
#Compute exact optimum allocation for PEVS-ADM 2016 data set for n = 50k
#Assumes tau = 1 by default
pevs_optE_alloc_50k <- opt_nh_nonresp(N_h = pevs_adm_2016_rrs$Nhat_h,
                                      phibar_h = pevs_adm_2016_rrs$rr_h,
                                      n_total = 50000)

#Merge results into data frame
(pevs_adm_alloc_50k_merged <- pevs_adm_2016_rrs \%>\%
  dplyr::mutate(n_h_optE = c(pevs_optE_alloc_50k),
                zeta_h_optE = attr(pevs_optE_alloc_50k,"zeta_h_prev")))



#For comparison purposes, compute approximate version of proposed allocation
pevs_optA_alloc_50k <- opt_nh_nonresp_oneiter(N_h = pevs_adm_2016_rrs$Nhat_h,
                                              phibar_h = pevs_adm_2016_rrs$rr_h,
                                              n_total = 50000)

#Merge results to previous tibble and reorder columns to be adjacent
pevs_adm_alloc_50k_merged \%>\%
  dplyr::mutate(n_h_optA = c(pevs_optA_alloc_50k)) \%>\%
  dplyr::relocate(zeta_h_optE, .after = "n_h_optA")

}
\references{
Mendelson, J., & Elliott, M. R. (in press).
\emph{Journal of Survey Statistics and Methodology}.
}
\keyword{allocation}
