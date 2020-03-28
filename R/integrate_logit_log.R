#' internal function used for integration of densities which appears
#' to be much more stable from -Inf to +Inf in the logit space while
#' the density to be integrated recieves inputs from 0 to 1 such that
#' the inverse distribution function must be used. The integral solved
#' is int_x dmix(mix,x) integrand(x) where integrand must be given as
#' log and we integrate over the support of mix.
#'
#' integrate density in logit space and split by component such
#' that the quantile function of each component is used. This
#' ensures that the R implementation of the quantile function is
#' always used.
#'
#' @param log_integrand function to integrate over which must return the log(f)
#' @param mix density over which to integrate
#' @param Lplower logit of lower cumulative density
#' @param Lpupper logit of upper cumulative density
#'
#' @keywords internal
integrate_density_log <- function(log_integrand, mix, Lplower=-Inf, Lpupper=Inf) {
    .integrand_comp <- function(mix_comp) {
        function(l) {
            u  <- inv_logit(l)
            lp <- log(u)
            lnp <- log(inv_logit(-l))
            exp(lp + lnp + log_integrand(qmix(mix_comp, u)))
        }
    }

    Nc <- ncol(mix)

    ## integrate by component of mix separatley to increase precision
    if(Nc == 1) {
        return(.integrate(.integrand_comp(mix), Lplower, Lpupper))
    } else {
        return(sum(vapply(1:Nc, function(comp) {
                              .integrate(.integrand_comp(mix[[comp, rescale=TRUE]]), Lplower, Lpupper)
                          }, c(0.1)) * mix[1,]))
    }
}

.integrate <- function(integrand, lower, upper) {
    integrate_args_user <- getOption("RBesT.integrate_args", list())
    args <- modifyList(list(lower=lower, upper=upper,
                            rel.tol=.Machine$double.eps^0.25,
                            abs.tol=.Machine$double.eps^0.25,
                            subdivisions=1000),
                       integrate_args_user)

    integrate(integrand,
              lower=args$lower,
              upper=args$upper,
              rel.tol=args$rel.tol,
              abs.tol=args$abs.tol,
              subdivisions=args$subdivisions)$value
}
