#' Decision Function for 1 Sample Operating Characteristics
#'
#' The function sets up a 1 sample one-sided decision function with an
#' arbitrary number of conditions which have all to be met.
#'
#' @param pc vector of critical cumulative probabilities of the
#' difference distribution.
#' @param qc vector of respective critical values of the difference
#' distribution. Must match the length of \code{pc}.
#' @param lower.tail logical value selecting if the threshold is a
#' lower or upper bound.
#'
#' @details The function creates a one-sided decision function which
#' takes two arguments. The first argument is expected to be a mixture
#' (posterior) distribution. This distribution is tested whether it
#' fulfills all the required threshold conditions specified with the
#' \code{pc} and \code{qc} arguments and returns 1 of all conditions
#' are met and 0 otherwise. Hence, for \code{lower.tail=TRUE}
#' condition \eqn{i} is equivalent to
#'
#' \deqn{P(x \leq q_{c,i}) > p_{c,i}}
#' 
#' and the decision function is implemented as indicator function on
#' the basis of the heavy-side step function \eqn{H} which is \eqn{0}
#' for \eqn{x \leq 0} and \eqn{1} for \eqn{x > 0}. As all conditions
#' must be met, the final indicator function returns
#' 
#' \deqn{\Pi_i H_i(P(x \leq q_{c,i}) - p_{c,i} ).}
#' 
#' When the second argument is set to \code{TRUE} a distance metric is
#' returned component wise as
#'
#' \deqn{ D_i = \log(P(p < q_{c,i})) - \log(p_{c,i}) .}
#'
#' These indicator functions can be used as input for 1-sample OC
#' calculations using \code{\link{oc1S}}.
#' 
#' @seealso oc1S
#'
#' @export
oc1Sdecision <- function(pc=0.975, qc=0, lower.tail=TRUE) {
    assert_that(length(pc) == length(qc))
    lpc <- log(pc)
    fun <- function(post, dist=FALSE) {
        test <- pmix(post, qc, lower.tail=lower.tail, log.p=TRUE) - lpc
        if(dist)
            return(test)
        as.numeric(all(test > 0))
    }
    attr(fun, "pc") <- pc
    attr(fun, "qc") <- qc
    attr(fun, "lower.tail") <- lower.tail
    class(fun) <- c("oc1Sdecision", "function")
    fun
}

#' @export
print.oc1Sdecision <- function(x, ...) {
    cat("1 sample decision function\n")
    cat("Conditions for acceptance:\n")
    qc <- attr(x, "qc")
    pc <- attr(x, "pc")
    low <- attr(x, "lower.tail")
    cmp <- ifelse(low, "<=", ">")
    for(i in seq_along(qc)) {
        cat(paste0("P(x ", cmp, " ", qc[i], ") > ", pc[i], "\n"))
    }
    invisible(x)
}
