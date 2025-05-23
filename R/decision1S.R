#' Decision Function for 1 Sample Designs
#'
#' The function sets up a 1 sample one-sided decision function with an
#' arbitrary number of conditions.
#'
#' @param pc Vector of critical cumulative probabilities.
#' @param qc Vector of respective critical values. Must match the length of `pc`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities
#' are \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#'
#' @details The function creates a one-sided decision function which
#' takes two arguments. The first argument is expected to be a mixture
#' (posterior) distribution. This distribution is tested whether it
#' fulfills all the required threshold conditions specified with the
#' `pc` and `qc` arguments and returns 1 if all conditions
#' are met and 0 otherwise. Hence, for `lower.tail=TRUE`
#' condition \eqn{i} is equivalent to
#'
#' \deqn{P(\theta \leq q_{c,i}) > p_{c,i}}
#'
#' and the decision function is implemented as indicator function on
#' the basis of the heavy-side step function \eqn{H(x)} which is \eqn{0}
#' for \eqn{x \leq 0} and \eqn{1} for \eqn{x > 0}. As all conditions
#' must be met, the final indicator function returns
#'
#' \deqn{\Pi_i H_i(P(\theta \leq q_{c,i}) - p_{c,i} ).}
#'
#' When the second argument is set to `TRUE` a distance metric is
#' returned component-wise per defined condition as
#'
#' \deqn{ D_i = \log(P(\theta < q_{c,i})) - \log(p_{c,i}) .}
#'
#' These indicator functions can be used as input for 1-sample
#' boundary, OC or PoS calculations using [oc1S()] or
#' [pos1S()] .
#'
#' @family design1S
#'
#' @return The function returns a decision function which takes two
#' arguments. The first argument is expected to be a mixture
#' (posterior) distribution which is tested if the specified
#' conditions are met. The logical second argument determines if the
#' function acts as an indicator function or if the function returns
#' the distance from the decision boundary for each condition in
#' log-space, i.e. the distance is 0 at the decision boundary,
#' negative for a 0 decision and positive for a 1 decision.
#'
#' @references Neuenschwander B, Rouyrre N, Hollaender H, Zuber E,
#' Branson M. A proof of concept phase II non-inferiority
#' criterion. *Stat. in Med.*. 2011, 30:1618-1627
#'
#' @examples
#'
#' # see Neuenschwander et al., 2011
#'
#' # example is for a time-to-event trial evaluating non-inferiority
#' # using a normal approximation for the log-hazard ratio
#'
#' # reference scale
#' s <- 2
#' theta_ni <- 0.4
#' theta_a <- 0
#' alpha <- 0.05
#' beta <- 0.2
#' za <- qnorm(1 - alpha)
#' zb <- qnorm(1 - beta)
#' n1 <- round((s * (za + zb) / (theta_ni - theta_a))^2) # n for which design was intended
#' nL <- 233
#' c1 <- theta_ni - za * s / sqrt(n1)
#'
#' # flat prior
#' flat_prior <- mixnorm(c(1, 0, 100), sigma = s)
#'
#' # standard NI design
#' decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
#'
#' # for double criterion with indecision point (mean estimate must be
#' # lower than this)
#' theta_c <- c1
#'
#' # double criterion design
#' # statistical significance (like NI design)
#' dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
#' # require mean to be at least as good as theta_c
#' dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
#' # combination
#' decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)
#'
#' theta_eval <- c(theta_a, theta_c, theta_ni)
#'
#' # we can display the decision function definition
#' decComb
#'
#' # and use it to decide if a given distribution fulfills all
#' # criterions defined
#' # for the prior
#' decComb(flat_prior)
#' # or for a possible outcome of the trial
#' # here with HR of 0.8 for 40 events
#' decComb(postmix(flat_prior, m = log(0.8), n = 40))
#'
#' @export
decision1S <- function(pc = 0.975, qc = 0, lower.tail = TRUE) {
  assert_that(length(pc) == length(qc))
  lpc <- log(pc)
  fun <- function(mix, dist = FALSE) {
    test <- pmix(mix, qc, lower.tail = lower.tail, log.p = TRUE) - lpc
    if (dist) {
      return(test)
    }
    as.numeric(all(test > 0))
  }
  attr(fun, "pc") <- pc
  attr(fun, "qc") <- qc
  attr(fun, "lower.tail") <- lower.tail
  class(fun) <- c("decision1S", "function")
  fun
}

#' @export
print.decision1S <- function(x, ...) {
  cat("1 sample decision function\n")
  cat("Conditions for acceptance:\n")
  qc <- attr(x, "qc")
  pc <- attr(x, "pc")
  low <- attr(x, "lower.tail")
  cmp <- ifelse(low, "<=", ">")
  for (i in seq_along(qc)) {
    cat(paste0("P(theta ", cmp, " ", qc[i], ") > ", pc[i], "\n"))
  }
  invisible(x)
}

#' @describeIn decision1S Deprecated old function name. Please use
#' `decision1S` instead.
#' @export
oc1Sdecision <- function(pc = 0.975, qc = 0, lower.tail = TRUE) {
  deprecated("oc1Sdecision", "decision1S")
  return(decision1S(pc, qc, lower.tail))
}
