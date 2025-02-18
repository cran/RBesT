% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decision2S_boundary.R
\name{decision2S_boundary}
\alias{decision2S_boundary}
\alias{decision2S_boundary.betaMix}
\alias{decision2S_boundary.normMix}
\alias{decision2S_boundary.gammaMix}
\title{Decision Boundary for 2 Sample Designs}
\usage{
decision2S_boundary(prior1, prior2, n1, n2, decision, ...)

\method{decision2S_boundary}{betaMix}(prior1, prior2, n1, n2, decision, eps, ...)

\method{decision2S_boundary}{normMix}(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-06,
  Ngrid = 10,
  ...
)

\method{decision2S_boundary}{gammaMix}(prior1, prior2, n1, n2, decision, eps = 1e-06, ...)
}
\arguments{
\item{prior1}{Prior for sample 1.}

\item{prior2}{Prior for sample 2.}

\item{n1, n2}{Sample size of the respective samples. Sample size \code{n1} must be greater than 0 while sample size \code{n2} must be greater or equal to 0.}

\item{decision}{Two-sample decision function to use; see \code{\link{decision2S}}.}

\item{...}{Optional arguments.}

\item{eps}{Support of random variables are determined as the
interval covering \code{1-eps} probability mass. Defaults to
\eqn{10^{-6}}.}

\item{sigma1}{The fixed reference scale of sample 1. If left
unspecified, the default reference scale of the prior 1 is assumed.}

\item{sigma2}{The fixed reference scale of sample 2. If left
unspecified, the default reference scale of the prior 2 is assumed.}

\item{Ngrid}{Determines density of discretization grid on which
decision function is evaluated (see below for more details).}
}
\value{
Returns a function with a single argument. This function
calculates in dependence of the outcome \eqn{y_2} in sample 2 the
critical value \eqn{y_{1,c}} for which the defined design will
change the decision from 0 to 1 (or vice versa, depending on the
decision function).
}
\description{
The \code{decision2S_boundary} function defines a 2 sample design
(priors, sample sizes, decision function) for the calculation of
the decision boundary. A function is returned which calculates the
critical value of the first sample \eqn{y_{1,c}} as a function of
the outcome in the second sample \eqn{y_2}. At the decision
boundary, the decision function will change between 0 (failure) and
1 (success) for the respective outcomes.
}
\details{
For a 2 sample design the specification of the priors, the
sample sizes and the decision function, \eqn{D(y_1,y_2)}, uniquely
defines the decision boundary

\deqn{D_1(y_2) = \max_{y_1}\{D(y_1,y_2) = 1\},}{D_1(y_2) = max_{y_1}{D(y_1,y_2) = 1},}

which is the critical value of \eqn{y_{1,c}} conditional on the
value of \eqn{y_2} whenever the decision \eqn{D(y_1,y_2)} function
changes its value from 0 to 1 for a decision function with
\code{lower.tail=TRUE} (otherwise the definition is \eqn{D_1(y_2) =
\max_{y_1}\{D(y_1,y_2) = 0\}}{D_1(y_2) = max_{y_1}{D(y_1,y_2) =
0}}). The decision function may change at most at a single critical
value for given \eqn{y_{2}} as only one-sided decision functions
are supported. Here, \eqn{y_2} is defined for binary and Poisson
endpoints as the sufficient statistic \eqn{y_2 = \sum_{i=1}^{n_2}
y_{2,i}} and for the normal case as the mean \eqn{\bar{y}_2 = 1/n_2
\sum_{i=1}^{n_2} y_{2,i}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{decision2S_boundary(betaMix)}: Applies for binomial model with a mixture
beta prior. The calculations use exact expressions.  If the
optional argument \code{eps} is defined, then an approximate method
is used which limits the search for the decision boundary to the
region of \code{1-eps} probability mass. This is useful for designs
with large sample sizes where an exact approach is very costly to
calculate.

\item \code{decision2S_boundary(normMix)}: Applies for the normal model with known
standard deviation \eqn{\sigma} and normal mixture priors for the
means. As a consequence from the assumption of a known standard
deviation, the calculation discards sampling uncertainty of the
second moment. The function has two extra arguments (with
defaults): \code{eps} (\eqn{10^{-6}}) and \code{Ngrid} (10). The
decision boundary is searched in the region of probability mass
\code{1-eps}, respectively for \eqn{y_1} and \eqn{y_2}. The
continuous decision function is evaluated at a discrete grid, which
is determined by a spacing with \eqn{\delta_2 =
\sigma_2/\sqrt{N_{grid}}}. Once the decision boundary is evaluated
at the discrete steps, a spline is used to inter-polate the
decision boundary at intermediate points.

\item \code{decision2S_boundary(gammaMix)}: Applies for the Poisson model with a gamma
mixture prior for the rate parameter.  The function
\code{decision2S_boundary} takes an extra argument \code{eps} (defaults to \eqn{10^{-6}}) which
determines the region of probability mass \code{1-eps} where the
boundary is searched for \eqn{y_1} and \eqn{y_2}, respectively.

}}
\examples{

# see ?decision2S for details of example
priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
# the success criteria is for delta which are larger than some
# threshold value which is why we set lower.tail=FALSE
successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
# the futility criterion acts in the opposite direction
futilityCrit <- decision2S(c(0.90), c(40), TRUE)

# success criterion boundary
successBoundary <- decision2S_boundary(priorP, priorT, 10, 20, successCrit)

# futility criterion boundary
futilityBoundary <- decision2S_boundary(priorP, priorT, 10, 20, futilityCrit)

curve(successBoundary(x), -25:25 - 49, xlab = "y2", ylab = "critical y1")
curve(futilityBoundary(x), lty = 2, add = TRUE)

# hence, for mean in sample 2 of 10, the critical value for y1 is
y1c <- futilityBoundary(-10)

# around the critical value the decision for futility changes
futilityCrit(postmix(priorP, m = y1c + 1E-3, n = 10), postmix(priorT, m = -10, n = 20))
futilityCrit(postmix(priorP, m = y1c - 1E-3, n = 10), postmix(priorT, m = -10, n = 20))

}
\seealso{
Other design2S: 
\code{\link{decision2S}()},
\code{\link{oc2S}()},
\code{\link{pos2S}()}
}
\concept{design2S}
