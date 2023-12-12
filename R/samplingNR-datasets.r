#' PEVS-ADM 2016 response rates by poststrata
#'
#'@description
#'Population and response rate information for the Federal Voting Assistance
#'Program's 2016 Post-Election Voting Survey
#'of U.S. Active Duty Military, based on secondary analysis of
#' public use data.
#'
#'@details
#'This table is from our paper's secondary analysis of survey data from
#'FVAP's 2016 Post-Election Voting Survey of U.S. Active Duty Military
#'(PEVS-ADM; FVAP, 2017a).
#'The sampling frame reflects U.S. Active Duty Military, as of
#'July 2016, based on administrative data from the
#'U.S. Defense Manpower Data Center (DMDC); see FVAP's technical
#'report (2017b) for details.
#'Our analysis exclusively examines the control group
#'from FVAP's response rate experiment (\eqn{n = 77{,}333});
#'we further omit \eqn{1{,}785} sample members who were identified as ineligible
#'via administative data prior to survey fielding; these
#'units were in the July 2016 sampling frame, but were identified via
#'September 2016 personnel records as no longer being on active duty
#'(e.g., due to separations or retirements).
#'Population sizes are estimated from the survey's
#'base weights for the remaining cases.
#'
#'Poststratum response propensities are base-weighted
#'proportions of invitees who provided an eligible, complete response to the
#'survey; the denominator includes all invited units
#'(as opposed to excluding ineligibles, as in AAPOR response rates),
#'so that resulting sample allocation calculations will reflect the
#'effects of sample loss on the costs of data collection.
#'Note that a different choice of denominator would not have much impact on
#'the estimated response rates, considering that an estimated
#'99.3\% of the invited population meets the study's
#'screener-based eligibility criteria.
#'
#' @format Tibble containing 91 observations and 8 variables. Includes:
#' \describe{
#'     \item{\code{h}}{(chr) Poststrata}
#'     \item{\code{service}}{(fct) Military service}
#'     \item{\code{paygrade}}{(fct) Pay grade category}
#'     \item{\code{age}}{(fct) Age category}
#'     \item{\code{region}}{(fct) Region}
#'     \item{\code{sex}}{(chr) Sex}
#'     \item{\code{Nhat_h}}{(dbl) Estimated population size in poststratum h}
#'     \item{\code{rr_h}}{(dbl) Base-weighted empirical response propensity in poststratum h}
#' }
#'@keywords data
#'@docType data
#'@references
#'Federal Voting Assistance Program. (2017a). 2015-2016 Active Duty Military dataset.
#'Federal Voting Assistance Program, U.S. Department of Defense.
#'\url{https://www.fvap.gov/uploads/FVAP/Surveys/PEV51601_Pub.dta} [Accessed 17 January 2023].
#'
#'Federal Voting Assistance Program (2017b). Post-Election Voting Surveys:
#'Active Duty Military: Technical report 2016.
#'Technical report, Federal Voting Assistance Program,
#'U.S. Department of Defense.
#'\url{https://www.fvap.gov/uploads/FVAP/Reports/PEVS_ADM_TechReport_Final.pdf}
#'[Accessed 17 January 2023].
#'@noMd
"pevs_adm_2016_rrs"

