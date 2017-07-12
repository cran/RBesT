#' Effective Sample Size for a Conjugate Prior
#' 
#' Calculates the Effective Sample Size (ESS) for a mixture prior. The
#' ESS indicates how many experimental units the prior is roughly
#' equivalent to.
#'
#' @param mix Prior (mixture of conjugate distributions).
#' @param method Selects the used method. Can be either \code{moment} (default) or \code{morita}.
# @param s For \code{morita} method large constant to ensure that the prior scaled by this value is vague (default 100); see Morita et al. (2008) for details.
#' @param ... Optional arguments applicable to specific methods.
#'
#' @details The ESS is calculated using the either a moments based
#' approach or the more sophisticated method by \emph{Morita et
#' al. (2008)}. The moments based method is the default method and
#' provides conservative estimates of the ESS.
#'
#' For the moments method the mean and standard deviation of the
#' mixture are calculated and then approximated by the conjugate
#' distribution with the same mean and standard deviation. For
#' conjugate distributions, the ESS is well defined. See the examples
#' for a step-wise calculation in the beta mixture case.
#'
#' The Morita method used here evaluates the mixture prior at the mode
#' instead of the mean as proposed originally by Morita. The method
#' may lead to very optimistic ESS values, especially if the mixture
#' contains many components.
#'
#' @template conjugate_pairs
#'
#' @references Morita S, Thall PF, Mueller P.
#' Determining the effective sample size of a parametric prior.
#' \emph{Biometrics} 2008;64(2):595-602.
#' 
#' @examples
#' # Conjugate Beta example
#' a <- 5
#' b <- 15
#' prior <- mixbeta(c(1, a, b))
#' 
#' ess(prior)
#' (a+b) 
#'
#' # Beta mixture example
#' bmix <- mixbeta(rob=c(0.2, 1, 1), inf=c(0.8, 10, 2))
#'
#' ess(bmix)
#' # moments method is equivalent to
#' # first calculate moments
#' bmix_sum <- summary(bmix)
#' # then calculate a and b of a matching beta
#' ab_matched <- ms2beta(bmix_sum["mean"], bmix_sum["sd"])
#' # finally take the sum of a and b which are equivalent
#' # to number of responders/non-responders respectivley
#' round(sum(ab_matched))
#' 
#' ess(bmix, method="morita")
#'
#' # Normal mixture example
#' nmix <- mixnorm(rob=c(0.5, 0, 2), inf=c(0.5, 3, 4), sigma=10)
#'
#' ess(nmix)
#'
#' ## the reference scale determines the ESS
#' sigma(nmix) <- 20
#' ess(nmix)
#' 
#' # Gamma mixture example
#' gmix <- mixgamma(rob=c(0.3, 20, 4), inf=c(0.7, 50, 10))
#'
#' ess(gmix) ## interpreted as appropriate for a Poisson likelihood (default)
#'
#' likelihood(gmix) <- "exp"
#' ess(gmix) ## interpreted as appropriate for an exponential likelihood
#' 
#'
#' @export
ess <- function(mix, method=c("moment", "morita"), ...) UseMethod("ess")
#' @export
ess.default <- function(mix, method=c("moment", "morita"), ...) "Unknown density"


calc_loc <- function(mix, loc=c("mode", "median", "mean")) {
    loc <- match.arg(loc)
    if(loc == "mode") {
        tol <- .Machine$double.eps^0.25
        locEst <- mixmode(mix)

        if(length(attr(locEst, "modes")) > 1)
            warning("Detected multiple modes.\nThe ESS is determined for the largest mode, but ESS concept is ill-defined for multi-modal distributions.")
    }
    if(loc == "median") {
        locEst <- qmix(mix, 0.5)
    }
    if(loc == "mean") {
        locEst <- summary(mix, NULL)["mean"]
    }
    names(locEst) <- NULL
    return(locEst)
}
## function to calculate mixture info of arbitrary density; input
## needed is the density function, gradient and hessian of the log
## density with respect to x (data)
mixInfo <- function(mix, x, dens, gradl, hessl) {
    p <- mix[1,]
    a <- mix[2,]
    b <- mix[3,]
    lp <- log(p)
    ldensComp <- dens(x, a, b, log=TRUE)
    ldensMix <- log_sum_exp(lp + ldensComp)
    lwdensComp <- lp + ldensComp - ldensMix
    dgl <- gradl(x,a,b)
    dhl <- (hessl(x,a,b) + dgl^2)
    ## attempt numerically more stable log calculations if possible,
    ## i.e. if all sings are the same
    if(all(dgl < 0) || all(dgl > 0)) {
        gsum <- exp(2*log_sum_exp(lwdensComp + log(abs(dgl))))
    } else {
        gsum <- (sum(exp(lwdensComp)*dgl))^2
    }
    if(all(dhl < 0) || all(dhl > 0)) {
        hsum <- sign(dhl[1]) * exp(log_sum_exp(lwdensComp + log(abs(dhl))))
    } else {
        hsum <- (sum(exp(lwdensComp)*dhl))
    }
    gsum - hsum
}

## function to calculate the gradient of the log mixture
## mixLogGrad <- function(mix, x, dens, gradl) {
##     p <- mix[1,]
##     a <- mix[2,]
##     b <- mix[3,]
##     densMix <- dmix(mix,x)
##     densComp <- dens(x, a, b)
##     dgl <- gradl(x,a,b)
##     sum(p*densComp*dgl)/densMix    
## }


# prior effective sample size ESS for Beta-mixture priors
# based on
# Morita, Thall, Mueller (MTM) 2008 Biometrics
# only difference: evaluated at mode of prior rather than at mean; and the flattened prior are derived with respect to the scale of 1 instead of being relative to the input scale
# SW: speedup by using analytical results and use of bisectioning search
#' @export
ess.betaMix <- function(mix, method=c("moment", "morita"), ..., s=100) {

    method <- match.arg(method)

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        res <- round(sum(ms2beta(smix["mean"], smix["sd"])))
        names(res) <- NULL
        return( res )
    }

    alphaP <- mix[2,]
    betaP <- mix[3,]

    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- betaMixInfo(mix, locEst)

    ESSmax <- ceiling(sum(alphaP+betaP)) * 2

    ## alpha and beta of "flattened" priors
    alphaP0 <- locEst / s
    betaP0 <- (1-locEst) / s

    ## MTM paper would follow this:
    ##priorN <- sum(mix[1,] * (alphaP + betaP))
    ##alphaP0 <- mode * priorN / s
    ##betaP0 <- (1-mode) * priorN / s

    ## we warn if any of the mixture components has a scale (n) which
    ## is less than 10/s such that the
    if(any(rowSums(mix[2:3,,drop=FALSE]) < 10/s )) {
        warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
    }

    ## expected 2nd derivative at mode
    ed2p <- function(m) {
        yn <- seq(0,m)
        ## negative 2nd log-derivative at mode
        info <- betaInfo(locEst, alphaP0 + yn, betaP0 + m - yn)
        ## prior predictive
        sum(info * dmix(preddist(mix,m), yn) )
    }

    ## function to search for change of sign
    ed2pDiff <- function(m) {
        deriv2.prior - ed2p(m)
    }

    return(unname(uniroot_int(ed2pDiff, c(0,ESSmax))))
}

## derivative of a single log-beta
betaLogGrad <- function(x,a,b) {
    lxm1 <- log1p(-x)
    lx <- log(x)
    - (b-1) * exp(- lxm1) + (a-1)*exp(- lx)
}

## second derivative of a single log-beta
betaLogHess <- function(x,a,b) {
    lxm1 <- log1p(-x)
    lx <- log(x)
    - (b-1) * exp(-2 * lxm1) - (a-1)*exp(-2 * lx)
}

betaMixInfo <- function(mix,x) {
    mixInfo(mix, x, dbeta, betaLogGrad, betaLogHess)
}

## info metric for a single beta, i.e. negative second derivative of log beta
betaInfo <- function(x,a,b) {
    -betaLogHess(x,a,b)
}


#' @export
ess.gammaMix <- function(mix, method=c("moment", "morita"), s=100, ...) {

    method <- match.arg(method)
    lik <- likelihood(mix)

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        coef <- ms2gamma(smix["mean"], smix["sd"])
        names(coef) <- NULL
        if(lik == "poisson")
            return(unname(round(coef[2])))
        if(lik == "exp")
            return(unname(round(coef[1])))
        stop("Unkown likelihood")
    }

    ## Morita method
    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- gammaMixInfo(mix, locEst)

    if(lik == "poisson") {
        meanPrior <- summary(mix)["mean"]
        names(meanPrior) <- NULL
        priorN <- mix[3,,drop=FALSE]

        ## function to search for change of sign
        ed2pDiff <- function(m) {
            deriv2.prior - ( gammaInfo(locEst, locEst/s + m * meanPrior, 1/s))
        }

        ESSmax <- ceiling(sum(mix[3,])) * 2
    }

    if(lik == "exp") {
        priorN <- mix[2,,drop=FALSE]
        ## function to search for change of sign
        ed2pDiff <- function(m) {
            deriv2.prior - gammaInfo(locEst, 1/s + m, 1/(s*locEst))
        }

        ESSmax <- ceiling(sum(mix[2,])) * 2
    }

    if(any(priorN < 10/s )) {
        warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
    }

    return(unname(uniroot_int(ed2pDiff, c(0,ESSmax))))
}

## respective functions for the gamma distribution
gammaLogGrad <- function(x,a,b) {
    (a-1)/x - b
}
gammaLogHess <- function(x,a,b) {
    -(a-1)/x^2
}
gammaInfo <- function(x,a,b) {
    -gammaLogHess(x,a,b)
}
gammaMixInfo <- function(mix,x) {
    mixInfo(mix, x, dgamma, gammaLogGrad, gammaLogHess)
}


#' @export
ess.normMix <- function(mix, method=c("moment", "morita"), sigma, s=100, ...) {

    method <- match.arg(method)

    if(missing(sigma)) {
        sigma <- RBesT::sigma(mix)
        message("Using default prior reference scale ", sigma)
    }
    assert_number(sigma, lower=0)
    tauSq <- sigma^2

    ## note: sigma is reassigned the sd's of each mixture component
    mu <- mix[2,]
    sigma <- mix[3,]
    sigmaSq <- sigma^2

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        res <- round(tauSq / smix["sd"]^2)
        return( unname(res) )
    }

    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- normMixInfo(mix, locEst)

    ESSmax <- ceiling(sum( (1-1/s) * tauSq/sigmaSq )) * 2

    ## "flattened" priors
    muP0 <- locEst
    sigmaP0Sq <- s * max(sigmaSq)

    ## difference of info at locEst and the expected 2nd derivative at
    ## locEst
    ed2pDiff <- function(m) {
        deriv2.prior - normInfo(locEst, muP0, sqrt(1/(m/tauSq + 1/sigmaP0Sq)))
    }

    return(unname(uniroot_int(ed2pDiff, c(0,ESSmax))))
}

## derivative of a single log-normal
normLogGrad <- function(x,mu,sigma) {
    -1 * (x-mu)/(sigma)^2
}

## second derivative of a single log-normal
normLogHess <- function(x,mu,sigma) {
    -1/sigma^2
}

normMixInfo <- function(mix,x) {
    mixInfo(mix, x, dnorm, normLogGrad, normLogHess)
}

## info metric for a single norm, i.e. negative second derivative of log norm
normInfo <- function(x,mean,sigma) {
    -normLogHess(x,mean,sigma)
}

