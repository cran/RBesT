#'
#' Transform `gMAP` to `draws` objects
#'
#' @description
#' Transform a `gMAP` object to a format supported by the
#' \pkg{posterior} package.
#'
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df
#' @aliases as_draws_rvars as_draws_list
#'
#' @param x A `gMAP` object.
#' @param variable A character vector providing the variables to
#'   extract.  By default, all variables are extracted.
#' @param regex Logical; Should variable be treated as a (vector of)
#'   regular expressions? Any variable in `x` matching at least
#'   one of the regular expressions will be selected. Defaults to
#'   `FALSE`.
#' @param inc_warmup Should warmup draws be included? Defaults to
#'   `FALSE`.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To subset iterations, chains, or draws, use the
#'   [posterior::subset_draws()] method after
#'   transforming the input object to a `draws` object.
#'
#' @seealso [posterior::draws()]
#'   [posterior::subset_draws()]
#'
#' @template example-start
#' @examples
#'
#' set.seed(34563)
#' map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
#'   family = binomial,
#'   data = AS,
#'   tau.dist = "HalfNormal", tau.prior = 1,
#'   beta.prior = 2
#' )
#'
#' post_AS <- as_draws(map_AS)
#'
#' @template example-stop
#'
#' @name draws-RBesT
NULL

#' @rdname draws-RBesT
#' @importFrom posterior as_draws
#' @method as_draws gMAP
#' @export
#' @export as_draws
as_draws.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix gMAP
#' @export
#' @export as_draws_matrix
as_draws_matrix.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_matrix,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_array
#' @method as_draws_array gMAP
#' @export
#' @export as_draws_array
as_draws_array.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_array,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_df
#' @method as_draws_df gMAP
#' @export
#' @export as_draws_df
as_draws_df.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_df,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_list
#' @method as_draws_list gMAP
#' @export
#' @export as_draws_list
as_draws_list.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars gMAP
#' @export
#' @export as_draws_rvars
as_draws_rvars.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_rvars,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    object_name = deparse(substitute(x)),
    ...
  )
}

#' @keywords internal
.as_draws_conversion <- function(
  x,
  draws_converter,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  object_name = deparse(substitute(x)),
  ...
) {
  stopifnot(inherits(x, "gMAP"))
  if (!is.logical(inc_warmup) || length(inc_warmup) != 1L || is.na(inc_warmup)) {
    stop("'inc_warmup' must be a single non-missing logical value.", call. = FALSE)
  }
  draws <- .gmap_get_stored_draws(x, inc_warmup = inc_warmup)
  if (is.null(draws)) {
    stop("The model does not contain posterior draws.", call. = FALSE)
  }
  # subset variables
  draws <- subset_draws(draws, variable = variable, regex = regex)
  if (posterior::ndraws(draws) == 0L) {
    object_name <- paste(object_name, collapse = " ")
    warning(
      "gMAP object \"",
      object_name,
      "\" does not contain any samples.",
      call. = FALSE
    )
  }
  draws_converter(draws)
}

#' @keywords internal
.gmap_get_stored_draws <- function(x, inc_warmup = FALSE) {
  if (is.null(x$draws)) {
    return(NULL)
  }

  if (!inc_warmup) {
    return(x$draws)
  }

  if (is.null(x$draws_warmup)) {
    stop(
      "Warmup draws were not stored. Refit with ",
      "options(RBesT.MC.save_warmup = TRUE) before calling ",
      "'inc_warmup = TRUE'."
      ,
      call. = FALSE
    )
  }

  posterior::as_draws_array(abind::abind(
    as.array(x$draws_warmup),
    as.array(x$draws),
    along = 1
  ))
}
