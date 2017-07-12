#' Operating Characteristics for 1 Sample Design
#'
#' Calculates the frequency at which the decision function is
#' evaluated to 1 under a specified true scenario and decision
#' criteria in a one-sample experiment.
#'
#' @param prior prior for analysis.
#' @param n sample size for the experiment.
#' @param decision one-sample decision function to use, see \code{\link{oc1Sdecision}}.
#' @param ... optional arguments.
#' @param eps support of random variables are determined as the
#' interval covering \code{1-eps} probability mass . Defaults to
#' \code{10^-6}.
#'
#' @details The operating characteristics calculate the frequency with
#' which the decision function is evaluated to 1 under the assumption
#' of a given true distribution of the data defined by
#' \eqn{\theta}. The specification of the prior, the sample size and
#' the decision function, \eqn{D(y)}, uniquely defines the
#' decision boundary
#' 
#' \deqn{y_c = \sup_y\{D(y) = 0\},}{y_c = sup_{y}{D(y) = 0},}
#'
#' which is the critical value whenever the decision \eqn{D(y)}
#' function changes its value from 0 to 1 for a decision function
#' with \code{lower.tail=TRUE} (otherwise the definition is \eqn{y_c =
#' \inf_{y}\{D(y) = 1\}}{y_c = inf_{y}{D(y) = 1}}). The decision
#' function may change at most at a single critical value as only
#' one-sided decision functions are supported. Here,
#' \eqn{y} is defined for binary and Poisson endpoints as the sufficient
#' statistic \eqn{y = \sum_{i=1}^{n} y_i} and for the normal
#' case as the mean \eqn{\bar{y} = 1/n \sum_{i=1}^n
#' y_i}.
#' 
#' Calling the \code{oc1S} function calculates the critical value
#' \eqn{y_c} and returns a function which can be used to query
#' the critical value or evaluate the desired frequency which is
#' evaluated as
#'
#' \deqn{ F(y_c|\theta). }
#'
#' \emph{Note:} Internally, the above integration is carried out
#' assuming that it suffices to assume a true distribution of the
#' mean.
#'
#' The convention for the critical value \eqn{y_c} depends on whether
#' a left (\code{lower.tail=TRUE}) or right-sided decision function
#' (\code{lower.tail=FALSE}) is used. For \code{lower.tail=TRUE} the
#' critical value \eqn{y_c} is the largest value for which the
#' decision is 1, \eqn{D(y \leq y_c) = 1}, while for
#' \code{lower.tail=FALSE} then \eqn{D(y > y_c) = 1} holds. This is
#' aligned with the cumulative density function definition within R
#' (see for example \code{\link{pbinom}}).
#'
#' @references Neuenschwander B, Rouyrre N, Hollaender H, Zuber E,
#' Branson M. A proof of concept phase II non-inferiority
#' criterion. \emph{Stat. in Med.}. 2011, 30:1618-1627
#' 
#' @return Returns a function which when called with one argument
#' \code{theta} will return the frequencies at which the decision
#' function is evaluated to 1. If called with no argument, then the
#' critical value \eqn{y_c} is returned. Note that the returned
#' function takes vectors arguments.
#'
#' @examples
#'
#' # see Neuenschwander et al., 2011
#' 
#' # example is for a time-to-event trial evaluationg non-inferiority
#' # using a normal approximation for the log-hazard ratio
#' 
#' # reference scale
#' s <- 2
#' theta_ni <- 0.4
#' theta_a <- 0
#' alpha <- 0.05
#' beta  <- 0.2
#' za <- qnorm(1-alpha)
#' zb <- qnorm(1-beta)
#' n1 <- round( (s * (za + zb)/(theta_ni - theta_a))^2 )  # n for which design was intended
#' nL <- 233
#' c1 <- theta_ni - za * s / sqrt(n1)
#'
#' # flat prior
#' prior <- mixnorm(c(1,0,100), sigma=s)
#' 
#' # standard NI design
#' decA <- oc1Sdecision(1 - alpha, theta_ni, lower.tail=TRUE)
#'
#' # for double criterion with indecision point (mean estimate must be
#' # lower than this)
#' theta_c <- c1
#' 
#' # double criterion design
#' # statistical significance (like NI design)
#' dec1 <- oc1Sdecision(1-alpha, theta_ni, lower.tail=TRUE)
#' # require mean to be at least as good as theta_c
#' dec2 <- oc1Sdecision(0.5, theta_c, lower.tail=TRUE)
#' # combination
#' decComb <- oc1Sdecision(c(1-alpha, 0.5), c(theta_ni, theta_c), lower.tail=TRUE)
#'
#' theta_eval  <- c(theta_a, theta_c, theta_ni)
#'
#' designA_n1  <- oc1S(prior, n1, decA)
#' designA_nL  <- oc1S(prior, nL, decA)
#' 
#' designC_n1  <- oc1S(prior, n1, decComb)
#' designC_nL  <- oc1S(prior, nL, decComb)
#' designC1_nL <- oc1S(prior, nL, dec1)
#' designC2_nL <- oc1S(prior, nL, dec2)
#'
#' designA_n1(theta_eval)
#' designA_nL(theta_eval)
#' designC_n1(theta_eval)
#' designC_nL(theta_eval)
#'
#' # note: the decision for the nL case is driven by the second
#' # criterion (for the mean) only as the critical value is lower
#' designC1_nL()
#' designC2_nL()
#' designC_nL()
#' 
#' @export
oc1S <- function(prior, n, decision, ...) UseMethod("oc1S")
#' @export
oc1S.default <- function(prior, n, decision, ...) "Unknown density"

#' @describeIn oc1S Applies for the beta-binomial model with a mixture
#' beta prior. The calculations use exact expressions.
#' @export
oc1S.betaMix <- function(prior, n, decision, ...) {

    VdecisionLazy <- Vectorize(function(r) { decision(postmix(prior, r=r, n=n)) - 0.25 } )

    ## find decision boundary
    bounds <- VdecisionLazy(c(0,n))
    lower.tail <- attr(decision, "lower.tail")
    if(prod(bounds) > 0) {
        ## decision is always the same
        if(lower.tail) {
            crit <- ifelse(bounds[1] < 0, 0, n+1)
        } else {
            crit <- ifelse(bounds[1] < 0, n, -1)
        }        
    } else {        
        crit <- uniroot_int(VdecisionLazy, c(0,n),
                            f.lower=bounds[1], f.upper=bounds[2])
    }

    ## crit is always pointing to the 0 just before the decision which
    ## is why we need a discrimination here
    if(lower.tail) {
        crit <- crit - 1
    }        
    
    design_fun <- function(theta) {
        if(missing(theta))
            return(crit)
        pbinom(crit, n, theta, lower.tail=lower.tail)
    }
    design_fun 
}


## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary1S_normMix <- function(decision, mix, n, lim) {

    sigma <- sigma(mix)
    
    cond_decisionStep <- function() {
        fn <- function(m) {
            decision(postmix(mix, m=m, se=sigma/sqrt(n))) - 0.75
        }
        Vectorize(fn)
    }

    uniroot(cond_decisionStep(), lim)$root
}

#' @describeIn oc1S Used for the normal-normal model with fixed
#' standard deviation (\eqn{\sigma}) of the sampling space and a
#' normal mixture prior. As a consequence from the mean assumption,
#' the calculation discards sampling uncertainty of the second
#' moment. The function has an extra argument \code{eps} (defaults to
#' \eqn{10^{-6}}). The critical value \eqn{y_c} is searched in the
#' region of probability mass \code{1-eps} for \eqn{y}.
#' @param sigma The fixed reference scale. If left unspecified, the
#' default reference scale of the prior is assumed.
#' @export
oc1S.normMix <- function(prior, n, decision, sigma, eps=1e-6, ...) {
    ## distributions of the means of the data generating distributions
    ## for now we assume that the underlying standard deviation
    ## matches the respective reference scales
    if(missing(sigma)) {
        sigma <- RBesT::sigma(prior)
        message("Using default prior reference scale ", sigma)
    }
    assert_number(sigma, lower=0)
    
    sd_samp <- sigma / sqrt(n)

    sigma(prior) <- sigma

    ## change the reference scale of the prior such that the prior
    ## represents the distribution of the respective means
    ##mean_prior <- prior
    ##sigma(mean_prior) <- sd_samp

    m <- summary(prior)["mean"]
    
    lim <- qnorm(p=c(eps/2, 1-eps/2), mean=m, sd=sd_samp)

    ## find the boundary of the decision function within the domain we integrate
    crit <- solve_boundary1S_normMix(decision, prior, n, lim)
    
    ## check where the decision is 1, i.e. left or right
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(theta) {
        if(missing(theta))
            return(crit)
        pnorm(crit, theta, sd_samp, lower.tail=lower.tail)
    }
    design_fun 
}


#' @describeIn oc1S Used for the Poisson-gamma case
#' (Poisson-exponential not yet supported). Takes an extra argument
#' \code{eps} (\eqn{10^6}) which determines the region of probability
#' mass \code{1-eps} where the boundary is searched for \eqn{y}.
#' @export
oc1S.gammaMix <- function(prior, n, decision, eps=1e-6, ...) {
    assert_that(likelihood(prior) == "poisson")

    cond_decisionStep <- function() {
        fn <- function(m) {
            decision(postmix(prior, n=n, m=m/n)) - 0.25
        }
        Vectorize(fn)
    }

    Vdecision <- cond_decisionStep()

    m <- summary(prior)["mean"]
    
    lambda_prior <- m * n
    
    lim <- qpois(p=c(eps/2, 1-eps/2), lambda=lambda_prior)
    
    bounds <- Vdecision(lim)
    ## make sure there is a decision somewhere
    while(prod(bounds) > 0) {
        lim[1] <- round(lim[1]/2)
        lim[2] <- round(lim[2]*2)
        bounds <- Vdecision(lim)
        if(diff(lim) > 1e9)
            break;
    }
    
    lower.tail <- attr(decision, "lower.tail")
    
    ## check if the decision is constantly 1 or 0
    if(prod(bounds) > 0) {
        ## decision is always the same
        if(lower.tail) {
            crit <- ifelse(bounds[1] < 0, 0, Inf)
        } else {
            crit <- ifelse(bounds[1] < 0, Inf, -1)
        }        
    } else {        
        ## find decision boundary
        crit <- uniroot_int(Vdecision,
                            c(lim[1],lim[2]),
                            f.lower=bounds[1],
                            f.upper=bounds[2])
    }

    ## crit is always pointing to the 0 just before the decision which
    ## is why we need a discrimination here
    if(lower.tail) {
        crit <- crit - 1
    }        

    ## crit is always pointing to the 0 just before the decision which
    ## is why we need a discrimination here
    design_fun <- function(theta) {
        if(missing(theta))
            return(crit)
        ppois(crit, n*theta, lower.tail=lower.tail)
    }
    design_fun 
}
