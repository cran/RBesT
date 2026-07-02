#' Resolve A Compact gMAP Fixture Spec Path
#'
#' Compact fixture specs live in `tests/testthat/fixtures-compact/` as R source
#' files. They are committed because they contain only recipe metadata and a
#' small parametric approximation, not full posterior draw objects.
#'
#' @param name Fixture name without the `_spec.R` suffix.
#' @return Absolute path to the compact fixture spec file.
compact_gmap_fixture_spec_path <- function(name) {
  testthat::test_path("fixtures-compact", paste0(name, "_spec.R"))
}

#' Load A Compact gMAP Test Fixture
#'
#' Compact fixtures build a draw-free `gMAP(..., chains = 0)` skeleton and then
#' inject posterior draws sampled from a stored MVN approximation. This is a
#' test-only pilot for lightweight pseudo-gMAP fixtures.
#'
#' @param name Fixture name without the `_spec.R` suffix.
#' @return A `gMAP` object with synthetic posterior draws.
load_compact_gmap_fixture <- function(name) {
  path <- compact_gmap_fixture_spec_path(name)
  if (!file.exists(path)) {
    stop("Missing compact gMAP fixture spec: ", path, call. = FALSE)
  }

  env <- new.env(parent = parent.frame())
  sys.source(path, envir = env)

  spec_name <- paste0(name, "_compact_spec")
  if (!exists(spec_name, envir = env, inherits = FALSE)) {
    stop("Compact gMAP fixture spec not found: ", spec_name, call. = FALSE)
  }

  inject_compact_gmap_draws(get(spec_name, envir = env, inherits = FALSE))
}

#' Fit A Compact MVN Mixture Spec From gMAP Draws
#'
#' This helper supports the pilot fixture workflow. It fits one MVN for the
#' model parameters and one MVN for the study-level theta draws.
#'
#' @param x A `gMAP` object with posterior draws.
#' @param seed Seed used when the compact fixture is rehydrated.
#' @param nc Number of MVN mixture components.
#' @return A compact fixture spec list.
create_compact_gmap_draw_spec <- function(x, seed, nc = 1L) {
  checkmate::assert_class(x, "gMAP")
  checkmate::assert_integerish(seed, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assert_integerish(nc, lower = 1, any.missing = FALSE, len = 1)

  theta_variables <- if (isTRUE(x$has_intercept)) {
    c("theta", "theta_pred")
  } else {
    "theta"
  }
  draws_theta <- posterior::as_draws_matrix(x, variable = theta_variables)
  draws_beta <- posterior::as_draws_matrix(x, variable = "beta")
  draws_tau <- posterior::as_draws_matrix(x, variable = "tau")
  draws_lp <- posterior::as_draws_matrix(x, variable = "lp__")
  tau_fixed <- compact_gmap_tau_is_fixed(draws_tau)
  tau_values <- if (tau_fixed) {
    stats::setNames(as.numeric(draws_tau[1, ]), colnames(draws_tau))
  } else {
    NULL
  }

  draws_model <- if (tau_fixed) {
    cbind(draws_beta, draws_lp)
  } else {
    draws_tau[draws_tau < .Machine$double.eps] <- .Machine$double.eps
    cbind(draws_beta, log(draws_tau), draws_lp)
  }

  mvn_model <- compact_gmap_mvn_mix(suppressWarnings(
    mixfit(draws_model, "mvnorm", Nc = nc)
  ))
  mvn_theta <- compact_gmap_mvn_mix(suppressWarnings(
    mixfit(draws_theta, "mvnorm", Nc = nc)
  ))

  list(
    name = NULL,
    seed = as.integer(seed),
    nc = as.integer(nc),
    has_intercept = isTRUE(x$has_intercept),
    tau_fixed = tau_fixed,
    tau_values = tau_values,
    variables_model = colnames(draws_model),
    variables_theta = colnames(draws_theta),
    variables = compact_gmap_draw_variables(
      draws_model,
      draws_theta,
      x,
      tau_variables = colnames(draws_tau)
    ),
    draws_dim = unname(dim(posterior::as_draws_array(x))),
    diag_variables = posterior::variables(x$draws_diag),
    metadata_mcmc = x$metadata_mcmc,
    mvn_model = mvn_model,
    mvn_theta = mvn_theta
  )
}

inject_compact_gmap_draws <- function(spec) {
  checkmate::assert_list(spec, names = "unique")
  checkmate::assert_function(spec$builder)

  skeleton <- spec$builder()
  if (!inherits(skeleton, "gMAP")) {
    stop("Compact fixture builder did not return a gMAP object.", call. = FALSE)
  }
  if (!identical(isTRUE(spec$has_intercept), isTRUE(skeleton$has_intercept))) {
    stop(
      "Compact gMAP fixture intercept setting does not match the skeleton.",
      call. = FALSE
    )
  }
  if (isTRUE(spec$tau_fixed)) {
    checkmate::assert_numeric(
      spec$tau_values,
      any.missing = FALSE,
      min.len = 1,
      names = "unique"
    )
  } else if (!is.null(spec$tau_values)) {
    stop(
      "Non-fixed tau compact gMAP fixtures must use `tau_values = NULL`.",
      call. = FALSE
    )
  }

  skeleton$draws <- rehydrate_compact_gmap_draws(spec, skeleton = skeleton)
  skeleton$draws_diag <- synthetic_compact_gmap_draws_diag(spec)
  skeleton$draws_warmup <- NULL
  skeleton$draws_warmup_diag <- NULL
  skeleton$metadata_mcmc <- spec$metadata_mcmc
  skeleton$thin <- 1L
  skeleton$backend <- "compact-fixture"

  beta_vars <- posterior::variables(
    posterior::as_draws_array(skeleton, variable = "beta")
  )
  tau_vars <- posterior::variables(
    posterior::as_draws_array(skeleton, variable = "tau")
  )

  beta_summary <- .gmap_summary(list(draws = skeleton$draws), variables = beta_vars)
  tau_summary <- .gmap_summary(list(draws = skeleton$draws), variables = tau_vars)
  skeleton$beta <- beta_summary[, "mean"]
  skeleton$tau <- tau_summary[, "mean"]
  names(skeleton$beta) <- colnames(skeleton$X)
  names(skeleton$tau) <- paste0("tau", seq_len(skeleton$n.tau.strata))

  sampler_summary <- .gmap_sampler_summary(
    list(draws = skeleton$draws),
    variables = c(beta_vars, tau_vars, "lp__")
  )
  skeleton$Rhat.max <- max(sampler_summary$rhat, na.rm = TRUE)

  attr(skeleton, "compact_fixture") <- list(
    name = spec$name,
    seed = spec$seed,
    nc = spec$nc,
    has_intercept = spec$has_intercept,
    tau_fixed = spec$tau_fixed,
    tau_values = spec$tau_values,
    variables_model = spec$variables_model
  )
  skeleton
}

rehydrate_compact_gmap_draws <- function(spec, skeleton) {
  checkmate::assert_class(skeleton, "gMAP")

  n_draws <- prod(spec$draws_dim[seq_len(2)])
  sampled <- withr::with_seed(
    spec$seed,
    list(
      model = rmix(spec$mvn_model, n_draws),
      theta = rmix(spec$mvn_theta, n_draws)
    )
  )

  model_draws <- sampled$model[, spec$variables_model, drop = FALSE]
  theta_draws <- sampled$theta[, spec$variables_theta, drop = FALSE]

  if (isTRUE(spec$tau_fixed)) {
    tau_draws <- matrix(
      rep(spec$tau_values, each = n_draws),
      nrow = n_draws,
      dimnames = list(NULL, names(spec$tau_values))
    )
  } else {
    tau_vars <- startsWith(colnames(model_draws), "tau[")
    tau_draws <- model_draws[, tau_vars, drop = FALSE]
    tau_draws[, ] <- exp(tau_draws)
  }

  predictive_draws <- NULL
  if (isTRUE(skeleton$has_intercept)) {
    theta_pred <- theta_draws[, "theta_pred"]
    theta_resp_pred <- skeleton$family$linkinv(theta_pred)
    theta_resp_pred <- compact_gmap_bound_response_draws(
      theta_resp_pred,
      skeleton$family
    )
    theta_pred <- matrix(theta_pred, ncol = 1, dimnames = list(NULL, "theta_pred"))
    theta_resp_pred <- matrix(
      theta_resp_pred,
      ncol = 1,
      dimnames = list(NULL, "theta_resp_pred")
    )
    predictive_draws <- cbind(
      theta_pred,
      theta_resp_pred
    )
  }

  draw_matrix <- cbind(
    theta_draws[, startsWith(colnames(theta_draws), "theta["), drop = FALSE],
    model_draws[, startsWith(colnames(model_draws), "beta["), drop = FALSE],
    tau_draws,
    predictive_draws,
    model_draws[, "lp__", drop = FALSE]
  )
  variables <- colnames(draw_matrix)
  dim(draw_matrix) <- spec$draws_dim
  dimnames(draw_matrix) <- list(
    iteration = seq_len(spec$draws_dim[[1]]),
    chain = seq_len(spec$draws_dim[[2]]),
    variable = variables
  )
  posterior::as_draws_array(draw_matrix)
}

synthetic_compact_gmap_draws_diag <- function(spec) {
  diag_array <- array(
    0,
    dim = c(spec$draws_dim[seq_len(2)], length(spec$diag_variables)),
    dimnames = list(
      iteration = seq_len(spec$draws_dim[[1]]),
      chain = seq_len(spec$draws_dim[[2]]),
      variable = spec$diag_variables
    )
  )

  diag_array[, , "accept_stat__"] <- 0.9
  diag_array[, , "stepsize__"] <- 0.01
  diag_array[, , "treedepth__"] <- 1
  diag_array[, , "n_leapfrog__"] <- 1
  diag_array[, , "divergent__"] <- 0
  diag_array[, , "energy__"] <- 1

  posterior::as_draws_array(diag_array)
}

compact_gmap_draw_variables <- function(draws_model, draws_theta, x, tau_variables) {
  theta_vars <- colnames(draws_theta)[startsWith(colnames(draws_theta), "theta[")]
  beta_vars <- colnames(draws_model)[startsWith(colnames(draws_model), "beta[")]
  c(
    theta_vars,
    beta_vars,
    tau_variables,
    if (x$has_intercept) c("theta_pred", "theta_resp_pred"),
    "lp__"
  )
}

compact_gmap_mvn_mix <- function(x) {
  attr(x, "df") <- NULL
  attr(x, "nobs") <- NULL
  attr(x, "lli") <- NULL
  attr(x, "Nc") <- NULL
  attr(x, "tol") <- NULL
  attr(x, "traceLli") <- NULL
  attr(x, "traceMix") <- NULL
  attr(x, "x") <- NULL
  class(x) <- c("mvnormMix", "mix")
  x
}

compact_gmap_tau_is_fixed <- function(draws_tau) {
  if (!ncol(draws_tau)) {
    return(FALSE)
  }
  all(apply(draws_tau, 2, function(x) {
    all(abs(x - x[[1]]) <= sqrt(.Machine$double.eps))
  }))
}

compact_gmap_bound_response_draws <- function(x, family) {
  eps <- .Machine$double.eps
  switch(
    family$link,
    logit = pmin(pmax(x, eps), 1 - eps),
    log = pmax(x, eps),
    x
  )
}
