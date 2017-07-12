
## EM for Beta Mixture Models (BMM) with Nc components

## TODO: revisit convergence criteria; consider changing to monitoring
## relative n and absolute mean; or revisit eps calculation with
## parameters mu and n.

EM_bmm_ab <- function(x, Nc, mix_init, verbose=TRUE, maxIter=500, tol, eps=0.1)
{
    N <- length(x)

    ## check data for 0 and 1 values which are problematic, but may be
    ## valid, depending on a and b. Moving these to eps or 1-eps
    ## ensures proper handling during fit.
    x0 <- x==0
    if(any(x0)) {
        message("Detected ", sum(x0), " value(s) which are exactly 0.\nTo avoid numerical issues during EM such values are moved to smallest eps on machine.")
        x[x0] <- .Machine$double.eps
    }
    x1 <- x==1
    if(any(x1)) {
        message("Detected ", sum(x1), " value(s) which are exactly 1.\nTo avoid numerical issues during EM such values are moved to one minus smallest eps on machine.")
        x[x1] <- 1-.Machine$double.eps
    }
    
    ## temporaries needed during EM
    Lx <- matrix(log(x), ncol=Nc, nrow=N)
    LxC <- matrix(log1p(-x), ncol=Nc, nrow=N)

    xRep <- rep(x, each=Nc)

    ## initialize randomly using KNN
    if(missing(mix_init)) {
        ##abmEst <- matrix(1+rlnorm(Nc*3, 0, log(5)/1.96), nrow=Nc)
        ##abmEst[,1] <- 1/Nc
        KNN <- suppressWarnings(knn(x, Nc, Niter.max=50))
        muInit <- as.vector(KNN$center)
        varInit <- tapply(x, KNN$cluster, var)
        nInit <- muInit*(1-muInit)/varInit - 1
        Nmax <- max(2, max(nInit))
        ## ensure n is positive for each cluster; if this is not the
        ## case, sample uniformly from the range of n we have
        ##Nneg <- nInit <= .Machine$double.eps
        Nsmall <- nInit <= 0.5
        if(any(Nsmall))
            nInit[Nsmall] <- runif(sum(Nsmall), 0.5, Nmax)
        nInitR <- 0.5 + rlnorm(Nc, log(nInit), log(5)/1.96)
        mixEst <- rbind(KNN$p, nInitR*muInit, nInitR*(1-muInit))
        rownames(mixEst) <- c("w", "a", "b")
    } else {
        mixEst <- mix_init
    }

    if(verbose) {
        cat("EM for beta mixture model.\n")
        cat("Initial estimates:\n")
        print(mixEst)
    }
    
    if(missing(tol)) {
        ## automatic selection of tolerance, based on targeted
        ## precision of alpha and beta (crude estimate)
        da <- digamma(mixEst[2,])
        db <- digamma(mixEst[3,])
        dab <- digamma(colSums(mixEst[2:3,,drop=FALSE]))
        tol <- min(abs((dab - da) * eps) + abs((dab - db) * eps))
        ##if(tol < 0.05)
        ##    warning("Tolerance set to very low value. Consider increasing eps, the precision of component estimates.")
    }
    
    if(verbose) {
        cat("Admissable change in log-likelihood tol =", tol,"\n")
    }

    lliLast <- -Inf
    iter <- 1
    logN <- log(N)
    traceMix <- list()
    traceLli <- c()

    ## find alpha and beta for a given component in log-space
    bmm_ml <- function(c1,c2) {
        function(par) {
            ab <- exp(par)
            s <- digamma(sum(ab))
            eq1 <- digamma(ab[1]) - s
            eq2 <- digamma(ab[2]) - s
            (eq1 - c1)^2 + (eq2 - c2)^2
        }
    }

    while(iter < maxIter) {
        ## calculate responsabilities from the likelihood terms;
        ## calculations are done in log-space to avoid numerical
        ## difficulties if some points are far away from some
        ## component and hence recieve very low density
        lli <- t(matrix(log(mixEst[1,]) + dbeta(xRep, mixEst[2,], mixEst[3,], log=TRUE), nrow=Nc))
        lnresp <- apply(lli, 1, log_sum_exp)
        ## the log-likelihood is then given by the sum of lresp norms
        lliCur <- sum(lnresp)
        lliDelta <- lliCur - lliLast
        ## record current state
        traceMix <- c(traceMix, list(mixEst))
        traceLli <- c(traceLli, lliCur)
        if(verbose) {
            cat("Iteration ", iter, ": log-likelihood = ", lliCur, "; delta = ", lliDelta, "\n", sep="")
        }
        if(iter != 1 & lliDelta < tol) {
            break
        }
        ## ... and the (log) responseability matrix follows from this by
        ## appropiate normalization.
        lresp <- sweep(lli, 1, lnresp, "-")
        resp <- exp(lresp)

        ## mean probability to be in a specific mixture component -> updates
        ## abmEst first colum
        lzSum <- apply(lresp, 2, log_sum_exp)
        zSum <- exp(lzSum)
        mixEst[1,] <- exp(lzSum - logN)

        c1 <- colSums(Lx * resp)/zSum
        c2 <- colSums(LxC * resp)/zSum

        ## now solve for new alpha and beta estimates jointly for each
        ## component
        for(i in 1:Nc) {
            theta <- c(log(mixEst[2:3,i]))
            Lest <- optim(theta, bmm_ml(c1[i], c2[i]))
            if(Lest$convergence != 0) {
                warning("Warning: Component", i, "in iteration", iter, "had convergence problems!")
            }
            ## TODO: consider adding the gradient to bmm_ml
            mixEst[2:3,i] <- exp(Lest$par)
        }

        lliLast <- lliCur
        iter <- iter + 1
    }
    if(iter == maxIter)
        warning("Maximum number of iterations reached.")

    mixEst <- mixEst[,order(mixEst[1,], decreasing=TRUE),drop=FALSE]
    colnames(mixEst) <- paste("comp", seq(Nc), sep="")
    class(mixEst) <- c("EM", "EMbmm", "betaMix")

    ## give further details
    attr(mixEst, "df") <- Nc-1 + 2*Nc
    attr(mixEst, "nobs") <- N
    attr(mixEst, "lli") <- lliCur

    attr(mixEst, "Nc") <- Nc
    
    attr(mixEst, "tol") <- tol
    attr(mixEst, "traceLli") <- traceLli
    attr(mixEst, "traceMix") <- lapply(traceMix, function(x) {class(x) <- "betaMix"; x})

    mixEst
}

