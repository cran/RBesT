#' Integrate a function against a mixture density
#'
#' Computes \eqn{\int g(x) f_{\text{mix}}(x) dx} where \eqn{f_{\text{mix}}}
#' is a mixture density. S3 generic dispatching on the mixture type.
#'
#' @param mix mixture density to integrate over (dispatch argument)
#' @param integrand function to integrate
#' @param ... additional arguments passed to methods
#'
#' @keywords internal
integrate_density <- function(mix, integrand, ...) UseMethod("integrate_density")

#' @exportS3Method
#' @describeIn integrate_density Default method using adaptive integration
#' @keywords internal
integrate_density.default <- function(
  mix,
  integrand,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6),
  ...
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) {
      u <- inv_logit(l)
      lp <- log_inv_logit(l)
      lnp <- log_inv_logit(-l)
      exp(lp + lnp) * integrand(qmix(mix_comp, u))
    }
  }
  Nc <- ncol(mix)

  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)

  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        ## ensure that the integrand is defined at the boundaries...
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        ## ... otherwise we avoid the boundaries by eps prob density:
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) integrand(x) * dmix(mix_comp, x),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

#' Integrate a log-space function against a mixture density
#'
#' Computes \eqn{\int \exp(\text{log\_integrand}(x)) f_{\text{mix}}(x) dx}.
#' Uses logit-space transformation for numerical stability with adaptive
#' integration. S3 generic dispatching on the mixture type.
#'
#' @param mix mixture density to integrate over (dispatch argument)
#' @param log_integrand function returning \code{log(g(x))}
#' @param ... additional arguments passed to methods
#'
#' @keywords internal
integrate_density_log <- function(mix, log_integrand, ...) UseMethod("integrate_density_log")

#' @exportS3Method
#' @describeIn integrate_density_log Default method using adaptive logit-space integration
#' @keywords internal
integrate_density_log.default <- function(
  mix,
  log_integrand,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6),
  ...
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) {
      u <- inv_logit(l)
      lp <- log_inv_logit(l)
      lnp <- log_inv_logit(-l)
      exp(lp + lnp + log_integrand(qmix(mix_comp, u)))
    }
  }

  Nc <- ncol(mix)

  ## integrate by component of mix separatley to increase precision
  ## when the density is not 0 at the boundaries integration, then
  ## the integration is performed on the natural scale. The check
  ## for that is done on the identity scale to avoid numerical
  ## issues.
  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)
  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) exp(log_integrand(x) + dmix(mix_comp, x, log = TRUE)),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

.integrate <- function(integrand, lower, upper) {
  integrate_args_user <- getOption("RBesT.integrate_args", list())
  args <- modifyList(
    list(
      lower = lower,
      upper = upper,
      rel.tol = .Machine$double.eps^0.25,
      abs.tol = .Machine$double.eps^0.25,
      subdivisions = 1000,
      stop.on.error = TRUE
    ),
    integrate_args_user
  )

  integrate(
    integrand,
    lower = args$lower,
    upper = args$upper,
    rel.tol = args$rel.tol,
    abs.tol = args$abs.tol,
    subdivisions = args$subdivisions,
    stop.on.error = args$stop.on.error
  )$value
}
