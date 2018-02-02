#' Beta-Binomial Probabilities
#' 
#' @param r,n number of successes (responders) out of n
#' @param a,b parameters of the Beta distribution for response probability
#' 
#' @details
#' r,n,a,b can be scalar or vectors. If vectors are used, they must be of the same length
#' 
#' @examples
#'
#' \dontrun{
#' # Ex1: Predictive distribution for uniform p
#' Ex1 =  dBetaBinomial( r=0:9,n=9,a=1,b=1)
#' Ex1
#' 
#' # Ex2: Predictive distribution at interim: n1=20, n=50
#' # Interim data: 4/20
#' # Probability to have 6 or more responders in 50 patients?
#' # That is: predictive probability >=2 in remaining 30?
#' 
#' # 1) Assume a weakly-informative Beta(a,1) prior with median 0.1 at trial start:
#' a = log(0.5)/log(0.1); b=1
#' p = dBetaBinomial(r=0:1,n=30,a=a+4,b=b+16)
#' 1-sum(p)
#' 
#' # 2) Assume a uniform prior at trial start:
#' p = dBetaBinomial(r=0:1,n=30,a=1+4,b=1+16)
#' 1-sum(p)
#' }
#' 
#' @keywords internal
#' 

## `dBetaBinomial` <- function(r,n,a,b,log=FALSE)
## {
##     assert_integerish(n, lower=0L, any.missing=FALSE)
##     assert_integerish(r, lower=0L, upper=max(n), any.missing=FALSE)
##     p <- lgamma(n+1)-lgamma(r+1)-lgamma(n-r+1)+lgamma(a+b)-lgamma(a)-lgamma(b)+lgamma(a+r)+lgamma(b+n-r)-lgamma(a+b+n);
##     if(log)
##         return(p)
##     exp(p);
## }


pBetaBinomial <- function(r, n, a, b, lower.tail=TRUE, log.p=FALSE) {
    return(.pBetaBinomial(r, n, a, b, lower.tail, log.p))
}