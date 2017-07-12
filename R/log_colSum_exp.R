#' Matrix version of log_sum_exp
#'
#' This version serves mostly for convenience as it is only slightly
#' faster then wrapping log_sum_exp into an apply call.
#'
#' @keywords internal
log_colSum_exp <- function(X) {
    nr <- nrow(X)
    if(nr==1) return(drop(X))
    nc <- ncol(X)
    Vxmax <- vector("numeric", 0)
    for(i in seq(nc)) {
        m <- which.max(X[,i])
        Vxmax[i] <- X[m,i]
        X[m,i] <- X[1,i]
    }
    log1p( .colSums(exp( X[-1,] - matrix(Vxmax, nr-1, nc, byrow=TRUE)  ), nr-1, nc) )+Vxmax
}
