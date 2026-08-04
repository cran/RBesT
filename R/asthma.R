#' Asthma exacerbation recurrent event data.
#'
#' Data set containing historical information for placebo arms of
#' relevant trials for the treatment of asthma. The primary outcome is
#' the rate of asthma exacerbations, a recurrent event modelled with a
#' negative binomial distribution. The full data set as published in
#' Holzhauer, Wang & Schmidli (2018) summarizes ten historical placebo
#' arms by their back-calculated log mean event rate and dispersion
#' together with the associated standard errors on the log scale.
#'
#' @format A data frame with 10 rows and 11 variables:
#' \describe{
#'   \item{study}{study label}
#'   \item{NCT}{ClinicalTrials.gov / ISRCTN registry identifier(s)}
#'   \item{d}{follow-up (exposure) duration in years}
#'   \item{n}{study size}
#'   \item{mu_hat}{estimated mean event rate}
#'   \item{log_mu_hat}{log mean event rate}
#'   \item{se_log_mu_hat}{standard error of the log mean event rate}
#'   \item{kappa_hat}{estimated dispersion parameter}
#'   \item{log_kappa_hat}{log dispersion parameter}
#'   \item{se_log_kappa_hat}{standard error of the log dispersion parameter}
#'   \item{phase}{development phase of the trial}
#' }
#'
#' @references \insertRef{holzhauer2018asthma}{RBesT}
#'
#' @template example-start
#' @examples
#' set.seed(34563)
#' asthma_ph3 <- subset(asthma, phase == "phase III")
#' map_asthma <- gMAP(cbind(log_mu_hat, se_log_mu_hat) ~ 1 + offset(log(d)) | study,
#'   family = gaussian,
#'   data = asthma_ph3,
#'   tau.dist = "HalfNormal", tau.prior = 0.5,
#'   beta.prior = cbind(0, 2)
#' )
#' @template example-stop
"asthma"
