#' Diagnostic plots for gMAP analyses
#'
#' @param x [gMAP()] object
#' @param size Controls size of forest plot.
#' @param linewidth Controls line sizes of traceplots.
#' @param ... Ignored.
#'
#' @details Creates MCMC diagnostics and a forest plot (including
#' model estimates) for a [gMAP()] analysis. For a
#' customized forest plot, please use the dedicated function
#' [forest_plot()].
#'
#' @template plot-help
#'
#' @return The function returns a list of [ggplot2::ggplot()]
#' objects.
#'
#' @method plot gMAP
#' @export
plot.gMAP <- function(x, size = NULL, linewidth = NULL, ...) {
  .gmap_warn_no_samples(x, object_name = deparse(substitute(x)))

  pl <- list()

  draws <- .gmap_draws_array(x)
  nuts_diag <- .gmap_diag_df(x$draws_diag)

  ## by default we return a small set of plots only
  plot_verbose <- getOption("RBesT.verbose", FALSE)

  div_opts <- list()
  if (!is.null(nuts_diag) && sum(nuts_diag$Value[nuts_diag$Parameter == "divergent__"]) > 0) {
    div_opts$np <- .gmap_all_diag_df(x)
    plot_verbose <- TRUE ## if any divergent transition happens,
    ## then we plot verbose in any case
  }

  tau_pars <- posterior::variables(.gmap_draws_array(x, variable = "tau"))
  beta_pars <- posterior::variables(.gmap_draws_array(x, variable = "beta"))

  tau_log_trans <- as.list(rep("log", length(tau_pars)))
  names(tau_log_trans) <- tau_pars

  if (plot_verbose) {
    ## traces are only shown if in verbose mode...
    draws_all <- .gmap_all_draws(x)
    n_warmup <- if (is.null(x$draws_warmup)) 0 else x$metadata_mcmc$warmup_saved
    pl$traceBeta <- do.call(
      bayesplot::mcmc_trace,
      c(
        list(
          x = draws_all,
          pars = beta_pars,
          size = size,
          n_warmup = n_warmup,
          facet_args = list(labeller = ggplot2::label_parsed)
        ),
        div_opts
      )
    ) +
      ggtitle(expression(paste("Trace of Regression Coefficient ", beta))) +
      bayesplot::facet_text(length(beta_pars) != 1) +
      xlab("Iteration")
    pl$traceTau <- do.call(
      bayesplot::mcmc_trace,
      c(
        list(
          x = draws_all,
          pars = tau_pars,
          size = size,
          n_warmup = n_warmup,
          facet_args = list(labeller = ggplot2::label_parsed)
        ),
        div_opts
      )
    ) +
      bayesplot::facet_text(length(tau_pars) != 1) +
      ggtitle(expression(paste("Trace of Heterogeneity Parameter ", tau))) +
      xlab("Iteration")
    pl$traceLogTau <- do.call(
      bayesplot::mcmc_trace,
      c(
        list(
          x = draws_all,
          pars = tau_pars,
          size = size,
          n_warmup = n_warmup,
          facet_args = list(labeller = ggplot2::label_parsed),
          transformations = tau_log_trans
        ),
        div_opts
      )
    ) +
      bayesplot::facet_text(length(tau_pars) != 1) +
      ggtitle(expression(paste(
        "Trace of Heterogeneity Parameter ",
        tau,
        " on log-scale"
      ))) +
      xlab("Iteration")

    ## ... as well as auxilary model parameters
    pl$densityBeta <- bayesplot::mcmc_dens_overlay(
      x = draws,
      pars = beta_pars,
      facet_args = list(
        labeller = ggplot2::label_parsed,
        strip.position = "bottom"
      )
    ) +
      ggtitle(expression(paste("Density of Regression Coefficient ", beta)))
    pl$densityTau <- bayesplot::mcmc_dens_overlay(
      x = draws,
      pars = tau_pars,
      facet_args = list(
        labeller = ggplot2::label_parsed,
        strip.position = "bottom"
      )
    ) +
      ggtitle(expression(paste("Density of Heterogeneity Parameter ", tau)))
    pl$densityLogTau <- bayesplot::mcmc_dens_overlay(
      x = draws,
      pars = tau_pars,
      facet_args = list(
        labeller = ggplot2::label_parsed,
        strip.position = "bottom"
      ),
      transformations = tau_log_trans
    ) +
      ggtitle(expression(paste(
        "Density of Heterogeneity Parameter ",
        tau,
        " on log-scale"
      )))
  }

  if (x$has_intercept) {
    pl$densityThetaStar <- bayesplot::mcmc_dens_overlay(
      x = draws,
      pars = "theta_resp_pred"
    ) +
      xlab(expression(theta[symbol("*")])) +
      bayesplot::facet_text(FALSE) +
      ggtitle(expression(paste("Density of MAP Prior ", theta[symbol("*")])))
    pl$densityThetaStarLink <- bayesplot::mcmc_dens_overlay(
      x = draws,
      pars = "theta_pred"
    ) +
      xlab(expression(theta[symbol("*")])) +
      bayesplot::facet_text(FALSE) +
      ggtitle(expression(paste(
        "Density of MAP Prior ",
        theta[symbol("*")],
        " (link scale)"
      )))

    pl$forest_model <- forest_plot(
      x,
      model = "both",
      size = if (is.null(size)) 1.25 else size,
      linewidth = if (is.null(linewidth)) 1.25 else linewidth
    )
  } else {
    message("No intercept defined.")
  }

  return(pl)
}
