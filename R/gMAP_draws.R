#' Extract non-permuted posterior draws from a Stan fit
#'
#' @keywords internal
.gmap_extract_draws_array <- function(fit, inc_warmup = FALSE) {
  posterior::as_draws_array(
    rstan::extract(
      fit,
      permuted = FALSE,
      inc_warmup = inc_warmup
    )
  )
}

#' Extract non-permuted sampler diagnostics from a Stan fit
#'
#' @keywords internal
.gmap_extract_diag_draws_array <- function(fit, inc_warmup = FALSE) {
  diag_list <- rstan::get_sampler_params(fit, inc_warmup = inc_warmup)

  if (!length(diag_list)) {
    return(NULL)
  }

  diag_vars <- colnames(diag_list[[1]])
  n_iter <- nrow(diag_list[[1]])
  n_chains <- length(diag_list)

  diag_array <- array(
    NA_real_,
    dim = c(n_iter, n_chains, length(diag_vars)),
    dimnames = list(
      iteration = as.character(seq_len(n_iter)),
      chain = as.character(seq_len(n_chains)),
      variable = diag_vars
    )
  )

  for (chain_id in seq_len(n_chains)) {
    diag_array[, chain_id, ] <- diag_list[[chain_id]][, diag_vars, drop = FALSE]
  }

  posterior::as_draws_array(diag_array)
}

#' Create An Empty Posterior Draws Array For A gMAP Skeleton
#'
#' @keywords internal
.gmap_empty_draws_array <- function(n_theta, n_beta, n_tau, has_intercept) {
  variables <- c(
    paste0("theta[", seq_len(n_theta), "]"),
    paste0("beta[", seq_len(n_beta), "]"),
    paste0("tau[", seq_len(n_tau), "]"),
    if (isTRUE(has_intercept)) c("theta_pred", "theta_resp_pred"),
    "lp__"
  )

  posterior::as_draws_array(array(
    numeric(),
    dim = c(0L, 0L, length(variables)),
    dimnames = list(
      iteration = character(),
      chain = character(),
      variable = variables
    )
  ))
}

#' Create An Empty Sampler Diagnostics Array For A gMAP Skeleton
#'
#' @keywords internal
.gmap_empty_diag_draws_array <- function() {
  variables <- c(
    "accept_stat__",
    "stepsize__",
    "treedepth__",
    "n_leapfrog__",
    "divergent__",
    "energy__"
  )

  posterior::as_draws_array(array(
    numeric(),
    dim = c(0L, 0L, length(variables)),
    dimnames = list(
      iteration = character(),
      chain = character(),
      variable = variables
    )
  ))
}

#' Create An NA Summary Table For Empty Draws
#'
#' @keywords internal
.gmap_empty_draws_summary <- function(variables, probs, row_names = NULL) {
  qnames <- names(posterior::quantile2(numeric(), probs = probs))
  out <- as.data.frame(matrix(
    NA_real_,
    nrow = length(variables),
    ncol = length(c("mean", "median", "sd", qnames)),
    dimnames = list(NULL, c("mean", "median", "sd", qnames))
  ))
  if (is.null(row_names)) {
    row_names <- variables
  }
  rownames(out) <- row_names
  out
}

#' Access Posterior Draws For gMAP Helpers
#'
#' @keywords internal
.gmap_draws_array <- function(x, variable = NULL, regex = FALSE) {
  if (inherits(x, "gMAP")) {
    draws <- .gmap_get_stored_draws(x)
    if (is.null(draws)) {
      stop("The model does not contain posterior draws.", call. = FALSE)
    }
    return(posterior::subset_draws(draws, variable = variable, regex = regex))
  }

  posterior::subset_draws(x$draws, variable = variable, regex = regex)
}

#' Thin Draws While Preserving Empty Skeleton Draws
#'
#' @keywords internal
.gmap_thin_draws <- function(x, thin) {
  checkmate::assert_integerish(
    thin,
    lower = 1,
    len = 1,
    any.missing = FALSE
  )
  thin <- as.integer(thin)

  if (posterior::ndraws(x) == 0L || thin <= 1L) {
    return(x)
  }

  posterior::thin_draws(x, thin = thin)
}

#' Warn When A gMAP Object Has No Posterior Samples
#'
#' @keywords internal
.gmap_warn_no_samples <- function(x, object_name) {
  draws <- .gmap_get_stored_draws(x)
  if (is.null(draws) || posterior::ndraws(draws) == 0L) {
    object_name <- paste(object_name, collapse = " ")
    warning(
      "gMAP object \"",
      object_name,
      "\" does not contain any samples.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' Collect MCMC metadata for a gMAP fit
#'
#' @keywords internal
.gmap_mcmc_metadata <- function(fit, thin_input, save_warmup) {
  warmup_saved <- fit@sim$warmup2
  if (length(warmup_saved)) {
    warmup_saved <- warmup_saved[[1]]
  } else {
    warmup_saved <- 0L
  }

  n_save <- fit@sim$n_save
  if (length(n_save)) {
    n_save <- n_save[[1]]
  } else {
    n_save <- 0L
  }

  list(
    iter = fit@sim$iter,
    warmup = fit@sim$warmup,
    warmup_saved = warmup_saved,
    chains = fit@sim$chains,
    n_save_per_chain = n_save,
    post_warmup_saved = n_save - warmup_saved,
    thin_input = thin_input,
    thin_post = 1L,
    save_warmup = isTRUE(save_warmup)
  )
}

#' Summarize selected posterior variables from stored draws
#'
#' @keywords internal
.gmap_summary <- function(
  x,
  variables,
  probs = c(0.025, 0.5, 0.975),
  transform = NULL,
  row_names = NULL
) {
  summary_draws <- .gmap_draws_array(x, variable = variables)
  if (posterior::ndraws(summary_draws) == 0L) {
    return(.gmap_empty_draws_summary(
      posterior::variables(summary_draws),
      probs = probs,
      row_names = row_names
    ))
  }

  if (!is.null(transform)) {
    summary_draws <- posterior::as_draws_matrix(
      transform(unclass(posterior::as_draws_matrix(summary_draws)))
    )
  }

  out <- posterior::summarise_draws(
    summary_draws,
    "mean",
    "median",
    "sd",
    ~posterior::quantile2(.x, probs = probs)
  )

  out <- as.data.frame(out)

  if (is.null(row_names)) {
    row_names <- out$variable
  }
  rownames(out) <- row_names
  out$variable <- NULL

  out
}

#' Summarize sampler diagnostics from stored draws
#'
#' @keywords internal
.gmap_sampler_summary <- function(x, variables = NULL) {
  draws <- .gmap_draws_array(x, variable = variables)
  if (posterior::ndraws(draws) == 0L) {
    return(data.frame(
      variable = posterior::variables(draws),
      rhat = NA_real_,
      ess_bulk = NA_real_,
      ess_tail = NA_real_
    ))
  }

  out <- posterior::summarise_draws(
    draws,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )

  names(out) <- sub("^posterior::", "", names(out))
  out
}

#' Combine warmup and post-warmup posterior draws
#'
#' @keywords internal
.gmap_all_draws <- function(x) {
  if (is.null(x$draws_warmup)) {
    return(x$draws)
  }

  posterior::bind_draws(
    x$draws_warmup,
    x$draws,
    along = "iteration"
  )
}

#' Convert diagnostic draws to bayesplot's long data frame format
#'
#' @keywords internal
.gmap_diag_df <- function(draws_diag) {
  if (is.null(draws_diag)) {
    return(NULL)
  }

  diag_vars <- posterior::variables(draws_diag)
  do.call(rbind, lapply(diag_vars, function(var) {
    df <- as.data.frame(
      posterior::as_draws_df(
        posterior::subset_draws(draws_diag, variable = var, scalar = TRUE)
      )
    )

    data.frame(
      Chain = df$.chain,
      Iteration = df$.iteration,
      Parameter = var,
      Value = df[[var]]
    )
  }))
}

#' Combine warmup and post-warmup diagnostics and return bayesplot format
#'
#' @keywords internal
.gmap_all_diag_df <- function(x) {
  if (is.null(x$draws_warmup_diag)) {
    return(.gmap_diag_df(x$draws_diag))
  }

  .gmap_diag_df(
    posterior::bind_draws(
      x$draws_warmup_diag,
      x$draws_diag,
      along = "iteration"
    )
  )
}

#' Count divergent transitions from stored sampler diagnostics
#'
#' @keywords internal
.gmap_divergence_count <- function(x) {
  if (is.null(x$draws_diag)) {
    return(0)
  }

  div <- posterior::subset_draws(x$draws_diag, variable = "divergent__")
  sum(as.array(div))
}
