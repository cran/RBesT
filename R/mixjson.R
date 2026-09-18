#' Write and Read a Mixture Object with JSON
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' These functions write and read a mixture object in the JSON format.
#'
#' @param mix A mixture object to be saved to JSON.
#' @param con A connection specifying where the JSON will be written
#'   to or read.
#' @param ... Additional arguments passed to the
#'   [jsonlite::toJSON()] and
#'   [jsonlite::fromJSON()] function for writing and
#'   reading, respectively.
#'
#' @details The mixture objects are written or read from the connection
#'   `con`, which can be a character string specifying a file
#'   path or a connection object as detailed in
#'   [base::connections()].
#'
#' When writing mixture objects as JSON it is strongly recommended to
#'   explicitly set the number of digits (argument `digits`) to be
#'   used for the numerical representation in order to control the
#'   accuracy of the JSON representation of the mixture object. If the
#'   mixture object inherits from the `"EM"` class (as is the
#'   case when the mixture is created using the
#'   [RBesT::mixfit()] function), then the mixture object
#'   will be cast to a simple mixture object such that diagnostics
#'   from the `"EM"` fitting procedure are dropped from the
#'   object. For easier readability the user is encouraged to set the
#'   argument `pretty=TRUE`, which is passed to the
#'   [jsonlite::toJSON()] function and makes the output more
#'   human readable.
#'
#' Note that when reading in mixture objects, then these are not
#' necessarily equal to the mixtures passed to the `write_mix_json`
#' function. This is a consequence of the limited precision of the
#' textual representation as defined by the `digits` argument.
#'
#' Whenever the requested precision is so low that a mixture weight is
#'   written as zero, then `write_mix_json` issues a warning, but
#'   still writes the mixture. Such a mixture can be read back
#'   without problems, since `read_mix_json` rescales the weights
#'   by default and zero weighted components are maintained. However,
#'   in case all mixture weights are written as zero, then
#'   `write_mix_json` aborts with an error as the mixture cannot be
#'   recovered any more from such a representation.
#'
#' @return The `write_mix_json` function does not return a value while
#'   the `read_mix_json` returns the mixture object stored in the
#'   connection specified.
#'
#' @family mixdist
#'
#' @examples
#' \dontrun{
#' nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
#'
#' write_mix_json(nm, "normal_mixture.json", pretty=TRUE, digits=1)
#'
#' mix <- read_mix_json("normal_mixture.json")
#' }
#' @name mixjson
NULL

#' @rdname mixjson
#' @export
write_mix_json <- function(mix, con, ...) {
  ## strip off all attributes which are added by the EM outputs and
  ## recast the class into a mixture only object
  if (inherits(mix, "EM")) {
    message("Dropping EM information from mixture object before serialization.")
    mix <- mixcombine(mix, weight = 1, rescale = FALSE)
  }
  dot_args <- list(...)
  if (!("digits" %in% names(dot_args))) {
    warning(
      "JSON serialization by default restricts number of digits to ",
      formals(jsonlite::toJSON)$digits,
      ".\nIt is recommended to set this option explicitly."
    )
  }
  umix <- unclass(mix)
  amix <- attributes(mix)
  amix$link <- amix$link$name
  check_json_weights(umix, dot_args$digits)
  json <- jsonlite::toJSON(list(meta = amix, comp = umix), ...)
  writeLines(json, con, useBytes = TRUE)
}

## determines the mixture weights as they are written to JSON by
## jsonlite and warns (errors) whenever some (all) weights become zero
check_json_weights <- function(umix, digits) {
  if (is.null(digits)) {
    digits <- formals(jsonlite::toJSON)$digits
  }
  weights <- umix[1, ]
  ## round-trip the weights through jsonlite in order to obtain the
  ## weights exactly as these are represented in the JSON output
  written <- jsonlite::fromJSON(jsonlite::toJSON(
    unname(weights),
    digits = digits
  ))
  is_zero <- written == 0
  if (!any(is_zero)) {
    return(invisible(TRUE))
  }
  if (all(is_zero)) {
    stop(
      "All mixture weights are written as zero with the requested precision (digits = ",
      digits,
      ").\nPlease increase the precision used for the JSON representation."
    )
  }
  labels <- colnames(umix)
  if (is.null(labels)) {
    labels <- which(is_zero)
  } else {
    labels <- labels[is_zero]
  }
  warning(
    "Mixture weight of component(s) ",
    paste0("'", labels, "'", collapse = ", "),
    " written as zero with the requested precision (digits = ",
    digits,
    ").\nThe stored mixture differs from the mixture given."
  )
  invisible(TRUE)
}

#' @rdname mixjson
#' @param rescale A logical value indicating whether to rescale the
#'   mixture weights so that they sum to 1. Defaults to `TRUE`.
#' @export
read_mix_json <- function(con, ..., rescale = TRUE) {
  json <- readLines(con)
  dot_args <- list(...)
  dot_args$simplifyMatrix <- FALSE
  mix_list <- do.call(jsonlite::fromJSON, c(list(txt = json), dot_args))
  mix <- do.call(
    rbind,
    lapply(mix_list$comp, function(component) unlist(component, use.names = FALSE))
  )
  meta <- mix_list$meta
  attributes(mix) <- meta
  attr(mix, "link") <- link_map[[mix_list$meta$link]]
  if (all(mix[1, ] == 0)) {
    stop(
      "All mixture weights are zero such that the mixture is not defined.\n",
      "Please increase the precision used when writing the mixture as JSON."
    )
  }
  if (rescale) {
    mix[1, ] <- mix[1, ] / sum(mix[1, ])
  }
  mix
}
