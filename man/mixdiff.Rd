% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixdiff.R
\name{mixdiff}
\alias{mixdiff}
\alias{dmixdiff}
\alias{pmixdiff}
\alias{qmixdiff}
\alias{rmixdiff}
\title{Difference of mixture distributions}
\usage{
dmixdiff(mix1, mix2, x)

pmixdiff(mix1, mix2, q, lower.tail = TRUE)

qmixdiff(mix1, mix2, p, lower.tail = TRUE)

rmixdiff(mix1, mix2, n)
}
\arguments{
\item{mix1}{first mixture density}

\item{mix2}{second mixture density}

\item{x}{vector of values for which density values are computed}

\item{q}{vector of quantiles for which cumulative probabilities are computed}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{P[X <= x]}, otherwise \eqn{P[X > x]}.}

\item{p}{vector of cumulative probabilities for which quantiles are computed}

\item{n}{size of random sample}
}
\value{
Respective density, quantile, cumulative density or random
numbers.
}
\description{
Density, cumulative distribution function, quantile
function and random number generation for the difference of two mixture
distributions.
}
\details{
If \eqn{x_1 \sim f_1(x_1)}{x_1 ~ f_1(x_1)} and \eqn{x_2 \sim
f_2(x_2)}{x_2 ~ f_2(x)}, the density of the difference \eqn{d
\equiv x_1 - x_2}{d = x_1 - x_2} is given by

\deqn{f_d(d) = \int f_1(u) \, f_2(u - d) \, du.}{f_d(d) = \int f_1(u) f_2(u - d)  du.}

The cumulative distribution function equates to

\deqn{F_d(d) = \int f_1(u) \, (1-F_2(u-d)) \, du.}{F_d(d) = \int f_1(u) (1-F_2(u-d)) du.}

Both integrals are performed over the full support of the
densities and use the numerical integration function
\code{\link[=integrate]{integrate()}}.
}
\examples{

# 1. Difference between two beta distributions, i.e. Pr( mix1 - mix2 > 0)
mix1 <- mixbeta(c(1, 11, 4))
mix2 <- mixbeta(c(1, 8, 7))
pmixdiff(mix1, mix2, 0, FALSE)

# Interval probability, i.e. Pr( 0.3 > mix1 - mix2 > 0)
pmixdiff(mix1, mix2, 0.3) - pmixdiff(mix1, mix2, 0)

# 2. two distributions, one of them a mixture
m1 <- mixbeta(c(1, 30, 50))
m2 <- mixbeta(c(0.75, 20, 50), c(0.25, 1, 1))

# random sample of difference
set.seed(23434)
rM <- rmixdiff(m1, m2, 1E4)

# histogram of random numbers and exact density
hist(rM, prob = TRUE, new = TRUE, nclass = 40)
curve(dmixdiff(m1, m2, x), add = TRUE, n = 51)

# threshold probabilities for difference, at 0 and 0.2
pmixdiff(m1, m2, 0)
mean(rM < 0)
pmixdiff(m1, m2, 0.2)
mean(rM < 0.2)

# median of difference
mdn <- qmixdiff(m1, m2, 0.5)
mean(rM < mdn)

# 95\%-interval
qmixdiff(m1, m2, c(0.025, 0.975))
quantile(rM, c(0.025, 0.975))

}
