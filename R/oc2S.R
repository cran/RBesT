#' Operating Characteristics for 2 Sample Design
#'
#' Calculates the frequency at which the decision function is
#' evaluated to 1 under a specified true scenario and decision
#' criteria in a two-sample experiment.
#'
#' @param prior1 prior for sample 1.
#' @param prior2 prior for sample 2.
#' @param n1,n2 sample size of the respective samples.
#' @param decision two-sample decision function to use, see \code{\link{oc2Sdecision}}.
#' @param ... optional arguments.
#' @param eps support of random variables are determined as the
#' interval covering \code{1-eps} probability mass . Defaults to
#' \code{10^-6}.
#' @param Ngrid determines density of discretization grid on which
#' decision function is evaluated (see below for more details).
#'
#' @details The operating characteristics calculate the frequency with
#' which the decision function is evaluated to 1 under the assumption
#' of a given true distribution of the data defined by \eqn{\theta_1}
#' and \eqn{\theta_2}. The specification of the priors, the sample
#' sizes and the decision function, \eqn{D(y_1,y_2)}, uniquely defines
#' the decision boundary
#' 
#' \deqn{D_1(y_2) = \sup_{y_1}\{D(y_1,y_2) = 1\},}{D_1(y_2) = sup_{y_1}{D(y_1,y_2) = 1},}
#'
#' which is the critical value of \eqn{y_{1,c}} conditional on the
#' value of \eqn{y_2} whenever the decision \eqn{D(y_1,y_2)} function
#' changes its value from 0 to 1 for a decision function with
#' \code{lower.tail=TRUE} (otherwise the definition is \eqn{D_1(y_2) =
#' \inf_{y_1}\{D(y_1,y_2) = 0\}}{D_1(y_2) = inf_{y_1}{D(y_1,y_2) =
#' 0}}). The decision function may change at most at a single critical
#' value as only one-sided decision functions are supported. Here,
#' \eqn{y_2} is defined for binary and Poisson endpoints as the
#' sufficient statistic \eqn{y_2 = \sum_{i=1}^{n_2} y_{2,i}} and for
#' the normal case as the mean \eqn{\bar{y}_2 = 1/n_2 \sum_{i=1}^{n_2}
#' y_{2,i}}.
#'
#' Calling the \code{oc2S} function calculates the decision boundary
#' \eqn{D_1(y_2)} and returns a function which can be used to evaluate
#' the decision boundary or evaluate the desired frequency which is
#' evaluated as
#'
#' \deqn{ \int f_2(y_2|\theta_2) F_1(D_1(y_2)|\theta_1) dy_2. }
#'
#' \emph{Note:} Internally, the above integration is carried out
#' assuming that it suffices to assume a true distribution of the
#' mean. 
#'
#' See below for examples and specifics for the supported mixture
#' priors.
#'
#' @return Returns a function which when called with two arguments
#' \code{theta1} and \code{theta2} will return the frequencies at
#' which the decision function is evaluated to 1. If called with a
#' named argument \code{y2} then the decision boundary \eqn{y_{1,c} =
#' D_1(y_2)} is returned. Note that the returned function takes
#' vector arguments.
#'
#' @references Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B.
#' Robust meta-analytic-predictive priors in clinical trials with historical control information.
#' \emph{Biometrics} 2014;70(4):1023-1032.
#'
#' @family oc2S
#'
#' @examples
#'
#' # example from Schmidli et al., 2015
#' dec <- oc2Sdecision(0.975, 0, lower.tail=FALSE)
#' 
#' N <- 40
#' prior_inf <- mixbeta(c(1, 4, 16))
#' prior_rob <- robustify(prior_inf, weight=0.2, mean=0.5)
#' prior_uni <- mixbeta(c(1, 1,  1))
#' N_ctl <- N - ess(prior_inf, method="morita")
#'
#' # compare designs with different priors
#' design_uni <- oc2S(prior_uni, prior_uni, N, N_ctl, dec)
#' design_inf <- oc2S(prior_uni, prior_inf, N, N_ctl, dec)
#' design_rob <- oc2S(prior_uni, prior_rob, N, N_ctl, dec)
#'
#' # decision boundary conditional on outcome of control
#' design_uni(y2=0:N_ctl)
#' design_inf(y2=0:N_ctl)
#' design_rob(y2=0:N_ctl)
#'
#' # type I error
#' curve(design_inf(x,x), 0, 1)
#' curve(design_uni(x,x), lty=2, add=TRUE)
#' curve(design_rob(x,x), lty=3, add=TRUE)
#' 
#' # power
#' curve(design_inf(0.2+x,0.2), 0, 0.5)
#' curve(design_uni(0.2+x,0.2), lty=2, add=TRUE)
#' curve(design_rob(0.2+x,0.2), lty=3, add=TRUE)
#' 
#'
#' @export
oc2S <- function(prior1, prior2, n1, n2, decision, ...) UseMethod("oc2S")
#' @export
oc2S.default <- function(prior1, prior2, n1, n2, decision, ...) "Unknown density"

#' @describeIn oc2S Applies for the beta-binomial model with a mixture
#' beta prior. The calculations use exact expressions. The decision
#' boundary is returned in alignment with the \code{lower.tail}
#' argument of the decision function, i.e. for \code{lower.tail=TRUE}
#' the returned critical value is the largest value \eqn{y_{1,c}} for
#' which the decision is 1 while for \code{lower.tail=FALSE} the
#' returned value is the largest value where the decision is still 0
#' which is compliant with the \code{\link{pbinom}} function. If the
#' extra argument \code{eps} is defined, then an approximate method is
#' used which limits the search for the decision boundary to the
#' region of \eqn{1-10^-6}. This is useful for designs with large
#' sample sizes where an exact approach is very costly to calculate.
#' @export
oc2S.betaMix <- function(prior1, prior2, n1, n2, decision, eps, ...) {
    ## only n2=0 is supported
    assert_number(n1, lower=1, finite=TRUE)
    assert_number(n2, lower=0, finite=TRUE)
    if(!missing(eps))
        assert_number(eps, lower=0, upper=0.1, finite=TRUE)

    if(!missing(eps) & (n1 > 1e9 | n2 > 1e9)) {
        warning("Large sample space. Consider setting eps=1e-6.")
    }

    cond_decisionDist <- function(post2cond) {
        fn <- function(m1) {
            ## Note: Subtracting from the decision 0.25 leads to
            ## negative decisions being at -0.25 while positives are
            ## at 0.75; since uniroot_int *always* returns the x which
            ## has lowest absolute value we are guaranteed that y2crit
            ## is just before the jump
            ##decision(post1cond, post2[[m2+1]]) - 0.25
            decision(postmix(prior1, r=m1, n=n1), post2cond) - 0.25
            ##decision(postmix(prior2, r=m2, n=n2), post1cond) - 0.25
        }
        Vectorize(fn)
    }

    ## saves the decision boundary conditional on the outcome of the
    ## second variable
    clim1 <- c(Inf, -Inf)
    clim2 <- c(Inf, -Inf)
    boundary <- c()
    full_boundary <- missing(eps)

    lower.tail <- attr(decision, "lower.tail")

    update_boundary <- function(lim1, lim2) {
        boundary <<- rep(NA, diff(lim2) + 1)
        clim2 <<- lim2
        clim1 <<- lim1
        for(y2 in lim2[1]:lim2[2]) {
            ## find decision point
            decFun <- cond_decisionDist(postmix(prior2, r=y2, n=n2))
            ind_llim <- decFun(lim1[1])
            ind_ulim <- decFun(lim1[2])
            y2ind <- y2 - lim2[1] + 1
            if(ind_llim < 0 & ind_ulim < 0) {
                ## then the decision is never 1
                boundary[y2ind] <<- -1
                next
            }
            if(ind_llim > 0 & ind_ulim > 0) {
                ## then the decision is always 1
                boundary[y2ind] <<- n1+1
                next
            }
            ## find boundary
            boundary[y2ind] <<- uniroot_int(decFun, lim1,
                                            f.lower=ind_llim,
                                            f.upper=ind_ulim)
        }
        if(lower.tail) {
            ## if lower.tail==TRUE, then the condition becomes true when
            ## going from large to small values, hence we need to integrate from
            ## 0 to boundary
            boundary <<- pmax(boundary - 1, -1)
        }
        boundary <- pmin(boundary, n1)
        return()
    }

    if(full_boundary)
        update_boundary(c(0, n1), c(0, n2))
    
    
    design_fun <- function(theta1, theta2, y2) {
        ## other-wise we calculate the frequencies at which the
        ## decision is 1 (probability mass with decision==1)

        ## in case n2==0, then theta2 is irrelevant
        if(n2 == 0 & missing(theta2))
            theta2 <- 0.5
        
        ## check if we need to recalculate the decision grid for the
        ## case of enabled approximate method
        if(!full_boundary) {
            if(!missing(y2)) {
                ## in case we want the critical values, setup theta1
                ## and theta2 such that y2 will always be in the
                ## domain in case they are not given
                if(missing(theta1))
                    theta1 <- c(0.1, 0.9)
                theta2 <- c(min(y2), max(y2))/n2
            }
            lim1 <- c(0,0)
            lim2 <- c(0,0)
            lim1[1] <- qbinom(  eps/2, n1, min(theta1))
            lim1[2] <- qbinom(1-eps/2, n1, max(theta1))
            lim2[1] <- qbinom(  eps/2, n2, min(theta2))
            lim2[2] <- qbinom(1-eps/2, n2, max(theta2))
            ## check if the decision grid needs to be recomputed
            if(lim1[1] < clim1[1] | lim1[2] > clim1[2] |
               lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
                ## ensure that lim1 never shrinks
                lim1[1] <- min(lim1[1], clim1[1])
                lim1[2] <- max(lim1[2], clim1[2])
                update_boundary(lim1, lim2)
            }
        }

        if(!missing(y2)) {
            ## make sure y2 is an integer which is the value of
            ## the second read-out for which we return the decision
            ## boundary
            ## TODO: handle case with eps with care
            assert_that(all(y2 == round(y2)))
            assert_that(all(y2 <= n2))
            assert_that(all(y2 >= 0))
            crit <- boundary[(y2 - clim2[1]) + 1]
            if(!full_boundary) {
                ## in case the lower boundary of the searched grid is not
                ## zero, then we cannot say anything about cases when the
                ## decision is always negative
                if(!lower.tail) {
                    ## in this case the decision changes from negative to
                    ## positive when going from small to large
                    ## values. Hence, if the decision is always negative,
                    ## then we can be sure of that we can never be sure,
                    ## but should the decision be negative at all values,
                    ## it can change at larger values.
                    crit[crit==n1+1] <- NA
                } else {
                    ## now the decision changes from positive to negative
                    ## when going from small to large => should the
                    ## decision not change in the clim1 domain then we do
                    ## not know if it happens later
                    if(clim1[1] > 0)
                        crit[crit==-1] <- NA
                    ## however, if crit==Inf then we can be sure that the
                    ## decision is indeed always positive
                }
            }
            return(crit)
        }
        

        T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names=NULL))
        if (inherits(T, "try-error")) {
            stop("theta1 and theta2 need to be of same size")
        }

        ## for each 0:n1 of the possible outcomes, calculate the
        ## probability mass past the boundary (log space) weighted with
        ## the density as the value for 1 occures (due to theta1)
        res <- matrix(-700, nrow=length(boundary), ncol=nrow(T))

        for(i in clim2[1]:clim2[2]) {
            y2ind <- i - clim2[1] + 1
            if(boundary[y2ind] == -1) {
                ## decision was always 0
                res[y2ind,] <- -700
            } else if(boundary[y2ind] == n1) {
                ## decision was always 1
                res[y2ind,] <- 0
            } else {
                ## calculate for all requested theta1 the probability mass
                ## past (or before) the boundary
                res[y2ind,] <- pbinom(boundary[y2ind], n1, T$theta1, lower.tail=lower.tail, log.p=TRUE)
            }
            ## finally weight with the density according to the occurence
            ## of i due to theta2; the pmax avoids -Inf in a case of Prob==0
            res[y2ind,] <- res[y2ind,] + pmax(dbinom(i, n2, T$theta2, log=TRUE), -700)
        }
        
        exp(log_colSum_exp(res))
    }
    design_fun
}


## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary2S_normMix <- function(decision, mix1, mix2, n1, n2, lim1, lim2, delta2) {
    grid <- seq(lim2[1], lim2[2], by=delta2)

    sigma1 <- sigma(mix1)
    sigma2 <- sigma(mix2)

    cond_decisionStep <- function(post2) {
        fn <- function(m1) {
            decision(postmix(mix1, m=m1, se=sigma1/sqrt(n1)), post2) - 0.75
        }
        Vectorize(fn)
    }
    
    Neval <- length(grid)
    ##cat("Calculating boundary from", lim2[1], "to", lim2[2], "with", Neval, "points\n")
    crit <- rep(NA, times=Neval)
    for(i in 1:Neval) {
        ind_fun <- cond_decisionStep(postmix(mix2, m=grid[i], se=sigma2/sqrt(n2)))
        dec_bounds <- ind_fun(lim1)
        ## if decision function is not different at boundaries, lim2
        ## is too narrow and we then enlarge
        while(prod(dec_bounds) > 0) {
            w <- diff(lim1)
            lim1 <- c(lim1[1] - w/2, lim1[2] + w/2)
            dec_bounds <- ind_fun(lim1)
        }
        crit[i] <- uniroot(ind_fun, lim1,
                           f.lower=dec_bounds[1], f.upper=dec_bounds[2])$root
    }
    
    cbind(grid, crit)
}

#' @describeIn oc2S Used for the normal-normal model with fixed
#' standard deviation (\eqn{\sigma}) of the sampling space and a
#' normal mixture prior. As a consequence from the mean assumption,
#' the calculation discards sampling uncertainty of the second
#' moment. The function has two extra arguments (with defaults):
#' \code{eps} (\eqn{10^{-6}}) and \code{Ngrid} (10). The decision
#' boundary is searched in the region of probability mass
#' \code{1-eps}, respectively for \eqn{y_1} and \eqn{y_2}. The
#' continuous decision function is evaluated at discrete steps which
#' are determined by \eqn{\delta_2 = \sigma_2/\sqrt{N_{grid}}}. Once the
#' decision boundary is evaluated at the discrete steps, a spline is
#' used to inter-polate the decision boundary at intermediate points.
#' @param sigma1 The fixed reference scale of sample 1. If left
#' unspecified, the default reference scale of the prior 1 is assumed.
#' @param sigma2 The fixed reference scale of sample 2. If left
#' unspecified, the default reference scale of the prior 2 is assumed.
#' @export
oc2S.normMix <- function(prior1, prior2, n1, n2, decision, sigma1, sigma2, eps=1e-6, Ngrid=10, ...) {
    ## distributions of the means of the data generating distributions
    ## for now we assume that the underlying standard deviation
    ## matches the respective reference scales
    if(missing(sigma1)) {
        sigma1 <- RBesT::sigma(prior1)
        message("Using default prior 1 reference scale ", sigma1)
    }
    if(missing(sigma2)) {
        sigma2 <- RBesT::sigma(prior2)
        message("Using default prior 2 reference scale ", sigma2)
    }
    assert_number(sigma1, lower=0)
    assert_number(sigma2, lower=0)
    
    sem1 <- sigma1 / sqrt(n1)
    sem2 <- sigma2 / sqrt(n2)

    sigma(prior1) <- sigma1
    sigma(prior2) <- sigma2

    ## only n2 can be zero
    assert_that(n1 >  0)
    assert_that(n2 >= 0)

    if(n2 == 0) sem2 <- sigma(prior2) / sqrt(1E-1)

    ## change the reference scale of the prior such that the prior
    ## represents the distribution of the respective means
    mean_prior1 <- prior1
    sigma(mean_prior1) <- sem1
    ##mean_prior2 <- prior2
    ##sigma(mean_prior2) <- sem2

    ## discretization step-size
    delta2 <- sem2/Ngrid

    ## the boundary function depends only on the samples sizes n1, n2,
    ## the priors and the decision, but not the assumed truths
    
    clim2 <- c(Inf, -Inf)

    ## the boundary function which gives conditional on the second
    ## variable the critical value where the decision changes
    boundary <- NA
    boundary_discrete <- matrix(NA, nrow=0, ncol=2)

    ## check where the decision is 1, i.e. left or right
    lower.tail <- attr(decision, "lower.tail")

    freq <- function(theta1, theta2) {
        lim2 <- qnorm(c(eps/2, 1-eps/2), theta2, sem2)
        if(n2 == 0) {
            return(pnorm(boundary(theta2), theta1, sem1, lower.tail=lower.tail))
            ##return(integrate(function(x) pnorm(boundary(x), theta1, sem1, lower.tail=lower.tail), lim2[1], lim2[2])$value)
        } else {
            
            return(integrate(function(x) dnorm(x, theta2, sem2) * pnorm(boundary(x), theta1, sem1, lower.tail=lower.tail), lim2[1], lim2[2], rel.tol=1E-10)$value)
        }
    }

    Vfreq <- Vectorize(freq)

    design_fun <- function(theta1, theta2, y2) {
        if(!missing(y2))
            theta2 <- y2

        ## in case n2==0, then theta2 is irrelevant
        if(n2 == 0 & missing(theta2))
            theta2 <- theta1
        
        lim2 <- c(qnorm(p=  eps/2, mean=min(theta2), sd=sem2)
                 ,qnorm(p=1-eps/2, mean=max(theta2), sd=sem2))

        ##cat("Boundaries", lim2, "cached", clim2,"\n")
        
        ## check if boundary function must be recomputed
        if(lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
            ## TODO: reuse old results and only compute what is extra needed
            ## find the boundary of the decision function within the domain we integrate
            ## note: the <<- assignment is needed to set the variable in the enclosure
            lim1 <- qmix(mean_prior1, c(eps/2, 1-eps/2))
            if(nrow(boundary_discrete) == 0) {
                ## boundary hasn't been calculated before, do it all
                boundary_discrete <<- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, lim2, delta2)
            } else {
                if(lim2[1] < clim2[1]) {
                    ## the lower bound is not low enough
                    boundary_extra <- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, c(lim2[1], max(lim2[1], clim2[1]-delta2)), delta2)
                    boundary_discrete <<- rbind(boundary_extra, boundary_discrete)
                }
                if(lim2[2] > clim2[2]) {
                    ## the upper bound is not large enough
                    boundary_extra <- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, c(min(lim2[2], clim2[2]+delta2), lim2[2]), delta2)
                    boundary_discrete <<- rbind(boundary_discrete, boundary_extra)
                }
            }
            boundary <<- splinefun(boundary_discrete[,1], boundary_discrete[,2])
            clim2 <<- lim2
        }

        if(!missing(y2)) {
            return(boundary(y2))
        }
        
        T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names=NULL))
        if (inherits(T, "try-error")) {
            stop("theta1 and theta2 need to be of same size")
        }

        do.call(Vfreq, T)
    }

    design_fun
}


#' @describeIn oc2S Used for the Poisson-gamma case
#' (Poisson-exponential not yet supported). Takes an extra argument
#' \code{eps} (\eqn{10^6}) which determines the region of probability
#' mass \code{1-eps} where the boundary is searched for \eqn{y_1} and
#' \eqn{y_2}, respectively.
#' @export
oc2S.gammaMix <- function(prior1, prior2, n1, n2, decision, eps=1e-6, ...) {
    assert_that(likelihood(prior1) == "poisson")
    assert_that(likelihood(prior2) == "poisson")

    # only the second n2 argument may be 0
    assert_that(n1 >  0)
    assert_that(n2 >= 0)
    
    cond_decisionStep <- function(post2) {
        fn <- function(m1) {
            decision(postmix(prior1, n=n1, m=m1/n1), post2) - 0.25
        }
        Vectorize(fn)
    }

    ## note: should be optimized by determining the critical point
    ## where the decision is taken conditional on the first only once

    clim1 <- c(Inf, -Inf)
    clim2 <- c(Inf, -Inf)
    boundary <- NA
    grid <- NA
    lower.tail <- attr(decision, "lower.tail")
    
    freq <- function(theta1, theta2) {
        lambda1 <- theta1 * n1
        lambda2 <- theta2 * n2

        exp(log_sum_exp(dpois(grid, lambda2, log=TRUE)
                        + ppois(boundary, lambda1, lower.tail=lower.tail, log.p=TRUE)))
    }

    Vfreq <- Vectorize(freq)

    design_fun <- function(theta1, theta2, y2) {

        ## in case n2==0, then theta2 is irrelevant
        if(n2 == 0 & missing(theta2))
            theta2 <- theta1

        if(!missing(y2)) {
            if(missing(theta1))
                theta1 <- summary(prior1)["mean"]
            lambda2 <- y2
        } else {
            lambda2 <- theta2 * n2
        }
        
        lambda1 <- theta1 * n1

        lim1 <- c(0, 0)
        lim1[1] <- qpois(  eps/2, min(lambda1))
        lim1[2] <- qpois(1-eps/2, max(lambda1))

        if(n2 == 0) {
            lim2 <- c(0,0)
        } else {
            lim2 <- c(0, 0)
            lim2[1] <- qpois(  eps/2, min(lambda2))
            lim2[2] <- qpois(1-eps/2, max(lambda2))
        }

        ## check if the boundary needs to be recomputed
        if(lim1[1] < clim1[1] | lim1[2] > clim1[2] |
           lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
            ## ensure that lim1 never shrinks
            lim1[1] <- min(lim1[1], clim1[1])
            lim1[2] <- max(lim1[2], clim1[2])
            grid <<- lim2[1]:lim2[2]
            Neval <- length(grid)
            boundary <<- rep(NA, Neval)
            for(i in 1:Neval) {
                cond_dec <- cond_decisionStep(postmix(prior2, n=n2, m=grid[i]/n2))
                low  <- cond_dec(lim1[1])
                high <- cond_dec(lim1[2])
                if(low < 0 & high < 0) {
                    boundary[i] <<- -1
                    next;
                }
                if(low > 0 & high > 0) {
                    boundary[i] <<- Inf
                    next;
                }
                boundary[i] <<- uniroot_int(cond_dec, lim1,
                                            f.lower=low,
                                            f.upper=high)
            }

            if(lower.tail) {
                ## if lower.tail==TRUE, then the condition becomes
                ## true when going from large to small values, hence
                ## we need to integrate from 0 to the boundary
                boundary <<- pmax(boundary - 1, -1)
            }

            ## save limits of new grid
            clim1 <<- lim1
            clim2 <<- lim2
        }

        if(!missing(y2)) {
            assert_that(all(y2 >= 0))
            crit <- boundary[y2 - clim2[1] + 1]
            ## in case the lower boundary of the searched grid is not
            ## zero, then we cannot say anything about cases when the
            ## decision is always negative
            if(!lower.tail) {
                ## in this case the decision changes from negative to
                ## positive when going from small to large
                ## values. Hence, if the decision is always negative,
                ## then we can be sure of that we can never be sure,
                ## but should the decision be negative at all values,
                ## it can change at larger values.
                crit[crit==Inf] <- NA
            } else {
                ## now the decision changes from positive to negative
                ## when going from small to large => should the
                ## decision not change in the clim1 domain then we do
                ## not know if it happens later
                if(clim1[1] > 0)
                    crit[crit==-1] <- NA
                ## however, if crit==Inf then we can be sure that the
                ## decision is indeed always positive
            }
            return(crit)
        }
        
        T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names=NULL))
        if (inherits(T, "try-error")) {
            stop("theta1 and theta2 need to be of same size")
        }
        do.call(Vfreq, T)
    }
    design_fun
}
