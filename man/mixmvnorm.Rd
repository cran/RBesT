% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixmvnorm.R
\name{mixmvnorm}
\alias{mixmvnorm}
\alias{msr2mvnorm}
\alias{print.mvnormMix}
\alias{summary.mvnormMix}
\alias{sigma.mvnormMix}
\title{Multivariate Normal Mixture Density}
\usage{
mixmvnorm(..., sigma, param = c("ms", "mn", "msr"))

msr2mvnorm(m, s, r, unlist = TRUE)

\method{print}{mvnormMix}(x, ...)

\method{summary}{mvnormMix}(object, ...)

\method{sigma}{mvnormMix}(object, ...)
}
\arguments{
\item{...}{List of mixture components.}

\item{sigma}{Reference covariance.}

\item{param}{Determines how the parameters in the list are
interpreted. See details.}

\item{m}{Mean vector.}

\item{s}{Standard deviation vector.}

\item{r}{Vector of correlations in column-major format of the lower
triangle of the correlation matrix.}

\item{unlist}{Logical. Controls whether the result is a flattened
vector (\code{TRUE}) or a list with mean \code{m} and covariance \code{s}
(\code{FALSE}). Defaults to \code{TRUE}.}

\item{object, x}{Multivariate normal mixture object.}
}
\value{
Returns a multivariate normal mixture with the specified
mixture components.
}
\description{
The multivariate normal mixture density and auxiliary
functions.
}
\details{
Each entry in the \code{...} argument list is a numeric
vector defining one component of the mixture multivariate
normal distribution. The first entry of the component defining
vector is the weight of the mixture component followed by the
vector of means in each dimension and finally a specification
of the covariance matrix, which depends on the chosen
parametrization. The covariance matrix is expected to be given
as numeric vector in a column-major format, which is standard
conversion applied to matrices by the vector concatenation
function \code{\link[base:c]{base::c()}}. Please refer to the examples
section below.

Each component defining vector can be specified in different ways
as determined by the \code{param} option:

\describe{
\item{ms}{Mean vector and covariance matrix \code{s}. Default.}
\item{mn}{Mean vector and number of observations. \code{n} determines
the covariance for each component via the relation \eqn{\Sigma/n}
with \eqn{\Sigma} being the known reference covariance.}
\item{msr}{Mean vector, standard deviations and correlations in
column-major format (corresponds to order when printing multi-variate
normal mixtures).}
}

The reference covariance \eqn{\Sigma} is the known covariance in
the normal-normal model (observation covariance). The function
\code{sigma} can be used to query the reference covariance and may
also be used to assign a new reference covariance, see examples
below. In case \code{sigma} is not specified, the user has to
supply \code{sigma} as argument to functions which require a
reference covariance.
}
\examples{

# default mean & covariance parametrization
S <- diag(c(1, 2)) \%*\% matrix(c(1, 0.5, 0.5, 1), 2, 2) \%*\% diag(c(1, 2))
mvnm1 <- mixmvnorm(
  rob = c(0.2, c(0, 0), diag(c(2, 2)^2)),
  inf = c(0.8, c(0.5, 1), S / 4), sigma = S
)

print(mvnm1)
summary(mvnm1)

set.seed(657846)
mixSamp1 <- rmix(mvnm1, 500)
colMeans(mixSamp1)

# alternative mean, sd and correlation parametrization
mvnm1_alt <- mixmvnorm(
  rob = c(0.2, c(0, 0), c(2, 2), 0.0),
  inf = c(0.8, c(0.5, 1), c(1, 2) / 2, 0.5),
  sigma = msr2mvnorm(s = c(1, 2), r = 0.5, unlist = FALSE)$s,
  param = "msr"
)

print(mvnm1_alt)

}
\seealso{
Other mixdist: 
\code{\link{mix}},
\code{\link{mixbeta}()},
\code{\link{mixcombine}()},
\code{\link{mixgamma}()},
\code{\link{mixjson}},
\code{\link{mixnorm}()},
\code{\link{mixplot}}
}
\concept{mixdist}
