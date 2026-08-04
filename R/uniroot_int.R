#' Find root of univariate function of integers
#'
#' Uses a bisectioning algorithm to search the give interval for a
#' change of sign and returns the integer which is closest to 0. When
#' `extendInt` is `"upX"` (f increasing through the root) or `"downX"`
#' (f decreasing), a bracket that does not enclose a sign change is
#' widened outward (doubling), bounded by `clamp`, until it brackets a
#' root. If no sign change exists within `clamp`, the constant sign is
#' reported as `-Inf` (f < 0) or `Inf` (f > 0); with `extendInt = "no"`
#' the legacy empty `numeric()` is returned instead.
#'
#' @keywords internal
uniroot_int <- function(
  f,
  interval,
  ...,
  f.lower = f(interval[1], ...),
  f.upper = f(interval[2], ...),
  extendInt = c("no", "upX", "downX"),
  clamp = interval,
  maxIter = 1000
) {
  extendInt <- match.arg(extendInt)
  assert_that(interval[1] < interval[2])

  lo <- interval[1]
  hi <- interval[2]
  fleft <- f.lower
  fright <- f.upper

  width <- max(hi - lo, 1)
  while (fleft * fright > 0 && extendInt != "no") {
    ## for "upX" the root is left of a positive bracket, right of a
    ## negative one; "downX" is the mirror image.
    grow_lo <- (extendInt == "upX") == (fleft > 0)
    if (grow_lo && lo > clamp[1]) {
      lo <- max(lo - width, clamp[1])
      fleft <- f(lo, ...)
    } else if (!grow_lo && hi < clamp[2]) {
      hi <- min(hi + width, clamp[2])
      fright <- f(hi, ...)
    } else {
      break
    }
    width <- width * 2
  }

  if (fleft * fright > 0) {
    ## constant sign across the (clamped) domain: no root. Report the
    ## sign so extend-mode callers can classify without re-evaluating f.
    if (extendInt == "no") {
      return(numeric())
    }
    return(if (fleft < 0) -Inf else Inf)
  }

  iter <- 0
  while ((hi - lo) > 1 & iter < maxIter) {
    mid <- floor((lo + hi) / 2)
    fmid <- f(mid, ...)
    if (fleft * fmid < 0) {
      hi <- mid
      fright <- fmid
    } else {
      lo <- mid
      fleft <- fmid
    }
    iter <- iter + 1
  }
  if (iter == maxIter) {
    warning("Maximum number of iterations reached.")
  }
  return(ifelse(abs(fleft) < abs(fright), lo, hi))
}

uniroot_int.all <- function(f, interval, maxIter = 1000, n = 100, ...) {
  assert_that(interval[1] < interval[2])

  xseq <- round(seq(interval[1], interval[2], len = n + 1))
  xseq <- xseq[!duplicated(xseq)]
  nu <- length(xseq) - 1
  mod <- f(xseq, ...)
  Equi <- xseq[which(mod == 0)]
  ss <- mod[1:nu] * mod[2:(nu + 1)]
  print(ss)
  ii <- which(ss < 0)
  print(ii)
  print(xseq[c(ii, ii[length(ii)] + 1)])
  for (i in ii) {
    Equi <- c(
      Equi,
      uniroot_int(f, c(xseq[i], xseq[i + 1]), ..., maxIter = maxIter)
    )
  }
  return(Equi)
}
