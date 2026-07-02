#' ---
#' title: "OC for Two-Stage Designs with Interim Analysis"
#' author: "Sebastian Weber"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' ---
#'
#' Demonstrates `oc2S_interim()`, a function for computing operating
#' characteristics of two-stage designs with an interim analysis (IA)
#' and a user-defined continuation rule (e.g., futility based on
#' probability of success).
#'
#' The approach avoids brute-force Monte Carlo by:
#' 1. Computing the **full-trial decision boundary once** (the expensive step).
#' 2. Using **Gauss-Hermite quadrature** over the interim data distribution
#'    to integrate out the IA outcomes deterministically.
#' 3. At each quadrature point, evaluating the conditional OC via the
#'    pre-computed boundary -- just `pnorm` + boundary lookup.
#'
#' This gives smooth, noise-free OC curves without nested MC.
#'
#' To render: `rmarkdown::render("oc2S_interim.R")`

library(RBesT)
library(MASS)
library(ggplot2)
theme_set(theme_bw())

## Load oc2S_interim() from inst/extra
source(system.file("extra", "oc2S_interim.R", package = "RBesT",
                    mustWork = TRUE))

#' ## Setup: same scenario as negbin_interim_futility vignette

## True parameters
lambda_ctrl <- 1.8
kappa_true <- 1.9
followup <- 1
log_exposure <- log(followup)

nb_family <- negative.binomial(theta = 1 / kappa_true)

## Per-patient sigma on log-rate scale
sigma_fn <- function(lambda, kappa, t) {
  sqrt((1 + kappa * lambda * t) / (lambda * t))
}
sigma_ctrl <- sigma_fn(lambda_ctrl, kappa_true, followup)

## Priors (non-informative)
uninf_ctrl <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)
uninf_treat <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)

## Sample sizes
n_treat <- 150
n_ctrl <- 150
n_treat_ia <- 45
n_ctrl_ia <- 45

## Final decision: P(delta < 0 | data) > 0.975
success_crit <- decision2S(0.975, 0, lower.tail = TRUE)

#' ## Define the IA rule
#'
#' Futility rule based on the posterior probability of a treatment
#' effect at the interim: continue if P(theta1 - theta2 < 0 | IA) > threshold.
#' This is cheap to evaluate (closed-form for normal mixtures) and
#' produces fast OC computation.

pos_threshold <- 0.10

ia_futility_rule <- function(post1_ia, post2_ia) {
  pmixdiff(post1_ia, post2_ia, 0) > pos_threshold
}

#' ## Compute OC curve

oc_ia <- oc2S_interim(
  prior1 = uninf_treat,
  prior2 = uninf_ctrl,
  n1 = n_treat,
  n2 = n_ctrl,
  n1_ia = n_treat_ia,
  n2_ia = n_ctrl_ia,
  decision = success_crit,
  ia_rule = ia_futility_rule,
  family = nb_family,
  offset1 = log_exposure,
  offset2 = log_exposure,
  Ngrid_ia = 21L
)

## Evaluate over a grid of treatment effects
log_mu_ctrl <- log(lambda_ctrl)
log_rr_grid <- seq(0, -0.8, by = -0.1)

cat("Computing OC curve...\n")
system.time(
  result <- oc_ia(
    theta1 = log_mu_ctrl + log_rr_grid,
    theta2 = rep(log_mu_ctrl, length(log_rr_grid))
  )
)

print(result)

#' ## Compare with fixed-design OC (no IA)

oc_fixed <- oc2S(
  uninf_treat, uninf_ctrl,
  n1 = n_treat, n2 = n_ctrl,
  decision = success_crit,
  family = nb_family,
  offset1 = log_exposure, offset2 = log_exposure
)

power_fixed <- oc_fixed(
  log_mu_ctrl + log_rr_grid,
  rep(log_mu_ctrl, length(log_rr_grid))
)

#' ## Plot

comp_df <- rbind(
  data.frame(log_rr = log_rr_grid, power = power_fixed,
             design = "No IA"),
  data.frame(log_rr = log_rr_grid, power = result$power,
             design = "IA (futility, grid)")
)

print(
  ggplot(comp_df, aes(log_rr, power, colour = design)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
    labs(
      x = "Log rate ratio (treatment vs control)",
      y = "Power",
      colour = NULL,
      title = "OC with futility IA (grid-based, no MC)"
    ) +
    scale_x_reverse() +
    theme(legend.position = "bottom")
)

#' ## Stopping probability

print(
  ggplot(result, aes(log_rr_grid, stop_prob)) +
    geom_line(linewidth = 0.8) +
    labs(
      x = "Log rate ratio",
      y = "P(stop for futility)",
      title = "Futility stopping probability (grid-based)"
    ) +
    scale_x_reverse() +
    coord_cartesian(ylim = c(0, 1))
)

#' ## Note on using pos2S-based rules
#'
#' For a PoS-based futility rule (more expensive but more
#' principled), replace the ia_rule with:
#' ```
#' ia_futility_pos <- function(post1_ia, post2_ia) {
#'   pos_fn <- pos2S(post1_ia, post2_ia,
#'                   n1 = n_treat - n_treat_ia,
#'                   n2 = n_ctrl - n_ctrl_ia,
#'                   decision = success_crit,
#'                   family = nb_family,
#'                   offset1 = log_exposure,
#'                   offset2 = log_exposure)
#'   pos_fn(post1_ia, post2_ia) >= 0.20
#' }
#' ```
#' This calls pos2S at each grid point and is slower (~minutes
#' with Ngrid_ia=21), but the outer integration remains exact.
#'
#' When using an informative prior only at the IA (e.g., MAP prior
#' for the control arm), pass `prior2_ia` to `oc2S_interim` and
#' use a 4-argument `ia_rule`:
#' ```
#' ia_futility_map <- function(post1_ia, post2_ia,
#'                             post1_ia_info, post2_ia_info) {
#'   # boundary uses uninformative posteriors (final analysis)
#'   pos_fn <- pos2S(post1_ia, post2_ia,
#'                   n1 = n_treat - n_treat_ia,
#'                   n2 = n_ctrl - n_ctrl_ia,
#'                   decision = success_crit,
#'                   family = nb_family,
#'                   offset1 = log_exposure,
#'                   offset2 = log_exposure)
#'   # prediction uses MAP-informed posterior
#'   pos_fn(post1_ia, post2_ia_info) >= 0.20
#' }
#'
#' oc_ia <- oc2S_interim(..., ia_rule = ia_futility_map,
#'                       prior2_ia = map_prior)
#' ```
