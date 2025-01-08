## ----SETTINGS-knitr, include=FALSE--------------------------------------------
## knitr settings used to build vignettes
library(RBesT)
library(knitr)
library(ggplot2)
theme_set(theme_bw())
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)
knitr::opts_chunk$set(
  dev = "ragg_png",
  dpi = 72,
  fig.retina = 2,
  fig.width = 1.62*4,
  fig.height = 4,
  fig.align = "center",
  out.width = "100%",
  pngquant = "--speed=1 --quality=50"
  )

## ----SETTINGS-sampling, include=FALSE-----------------------------------------
## sampling settings used to build vignettes
## setup up fast sampling when run on CRAN
is_CRAN <- Sys.getenv("NOT_CRAN", "true") != "true"
## NOTE: for running this vignette locally, please uncomment the
## following line:
## is_CRAN <- FALSE
.user_mc_options <- list()
if (is_CRAN) {
    .user_mc_options <- options(RBesT.MC.warmup=250, RBesT.MC.iter=500, RBesT.MC.chains=2, RBesT.MC.thin=1, RBesT.MC.control=list(adapt_delta=0.9))
}
set.seed(6475863)

## ----results="asis",echo=FALSE------------------------------------------------
kable(AS)

## -----------------------------------------------------------------------------
# load R packages
library(RBesT)
library(ggplot2)
theme_set(theme_bw()) # sets up plotting theme

set.seed(34563)
map_mcmc <- gMAP(cbind(r, n - r) ~ 1 | study,
  data = AS,
  tau.dist = "HalfNormal",
  tau.prior = 1,
  beta.prior = 2,
  family = binomial
)
print(map_mcmc)

## a graphical representation of model checks is available
pl <- plot(map_mcmc)

## a number of plots are immediately defined
names(pl)

## forest plot with model estimates
print(pl$forest_model)

## -----------------------------------------------------------------------------
set.seed(36546)
map_mcmc_sens <- update(map_mcmc, tau.prior = 1 / 2)
print(map_mcmc_sens)

## -----------------------------------------------------------------------------
map <- automixfit(map_mcmc)
print(map)
plot(map)$mix

## -----------------------------------------------------------------------------
round(ess(map, method = "elir")) ## default method
round(ess(map, method = "moment"))
round(ess(map, method = "morita"))

## -----------------------------------------------------------------------------
## add a 20% non-informative mixture component
map_robust <- robustify(map, weight = 0.2, mean = 1 / 2)
print(map_robust)

round(ess(map_robust))

## -----------------------------------------------------------------------------
ess_weight <- data.frame(weight = seq(0.05, 0.95, by = 0.05), ess = NA)
for (i in seq_along(ess_weight$weight)) {
  ess_weight$ess[i] <- ess(robustify(map, ess_weight$weight[i], 0.5))
}
ess_weight <- rbind(
  ess_weight,
  data.frame(
    weight = c(0, 1),
    ess = c(ess(map), ess(mixbeta(c(1, 1, 1))))
  )
)

ggplot(ess_weight, aes(weight, ess)) +
  geom_point() +
  geom_line() +
  ggtitle("ESS of robust MAP for varying weight of robust component") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 40, by = 5))

## -----------------------------------------------------------------------------
theta <- seq(0.1, 0.95, by = 0.01)
uniform_prior <- mixbeta(c(1, 1, 1))
treat_prior <- mixbeta(c(1, 0.5, 1)) # prior for treatment used in trial
lancet_prior <- mixbeta(c(1, 11, 32)) # prior for control   used in trial
decision <- decision2S(0.95, 0, lower.tail = FALSE)

design_uniform <- oc2S(uniform_prior, uniform_prior, 24, 6, decision)
design_classic <- oc2S(uniform_prior, uniform_prior, 24, 24, decision)
design_nonrobust <- oc2S(treat_prior, map, 24, 6, decision)
design_robust <- oc2S(treat_prior, map_robust, 24, 6, decision)

typeI_uniform <- design_uniform(theta, theta)
typeI_classic <- design_classic(theta, theta)
typeI_nonrobust <- design_nonrobust(theta, theta)
typeI_robust <- design_robust(theta, theta)

ocI <- rbind(
  data.frame(theta = theta, typeI = typeI_robust, prior = "robust"),
  data.frame(theta = theta, typeI = typeI_nonrobust, prior = "non-robust"),
  data.frame(theta = theta, typeI = typeI_uniform, prior = "uniform"),
  data.frame(theta = theta, typeI = typeI_classic, prior = "uniform 24:24")
)

ggplot(ocI, aes(theta, typeI, colour = prior)) +
  geom_line() +
  ggtitle("Type I Error")

## -----------------------------------------------------------------------------
summary(map)

## -----------------------------------------------------------------------------
ggplot(ocI, aes(theta, typeI, colour = prior)) +
  geom_line() +
  ggtitle("Type I Error - response rate restricted to plausible range") +
  coord_cartesian(xlim = c(0, 0.5))

## -----------------------------------------------------------------------------
delta <- seq(0, 0.7, by = 0.01)
mean_control <- summary(map)["mean"]
theta_active <- mean_control + delta
theta_control <- mean_control + 0 * delta

power_uniform <- design_uniform(theta_active, theta_control)
power_classic <- design_classic(theta_active, theta_control)
power_nonrobust <- design_nonrobust(theta_active, theta_control)
power_robust <- design_robust(theta_active, theta_control)

ocP <- rbind(
  data.frame(theta_active, theta_control, delta = delta, power = power_robust, prior = "robust"),
  data.frame(theta_active, theta_control, delta = delta, power = power_nonrobust, prior = "non-robust"),
  data.frame(theta_active, theta_control, delta = delta, power = power_uniform, prior = "uniform"),
  data.frame(theta_active, theta_control, delta = delta, power = power_classic, prior = "uniform 24:24")
)

ggplot(ocP, aes(delta, power, colour = prior)) +
  geom_line() +
  ggtitle("Power")

## -----------------------------------------------------------------------------
find_delta <- function(design, theta_control, target_power) {
  uniroot(
    function(delta) {
      design(theta_control + delta, theta_control) - target_power
    },
    interval = c(0, 1 - theta_control)
  )$root
}

target_effect <- data.frame(
  delta = c(
    find_delta(design_nonrobust, mean_control, 0.8),
    find_delta(design_classic, mean_control, 0.8),
    find_delta(design_robust, mean_control, 0.8),
    find_delta(design_uniform, mean_control, 0.8)
  ),
  prior = c("non-robust", "uniform 24:24", "robust", "uniform")
)

knitr::kable(target_effect, digits = 3)

## -----------------------------------------------------------------------------
## Critical values at which the decision flips are given conditional
## on the outcome of the second read-out; as we like to have this as a
## function of the treatment group outcome, we flip label 1 and 2
decision_flipped <- decision2S(0.95, 0, lower.tail = TRUE)
crit_uniform <- decision2S_boundary(uniform_prior, uniform_prior, 6, 24, decision_flipped)
crit_nonrobust <- decision2S_boundary(map, treat_prior, 6, 24, decision_flipped)
crit_robust <- decision2S_boundary(map_robust, treat_prior, 6, 24, decision_flipped)
treat_y2 <- 0:24
## Note that -1 is returned to indicated that the decision is never 1
ocC <- rbind(
  data.frame(y2 = treat_y2, y1_crit = crit_robust(treat_y2), prior = "robust"),
  data.frame(y2 = treat_y2, y1_crit = crit_nonrobust(treat_y2), prior = "non-robust"),
  data.frame(y2 = treat_y2, y1_crit = crit_uniform(treat_y2), prior = "uniform")
)

ggplot(ocC, aes(y2, y1_crit, colour = prior)) +
  geom_step() +
  ggtitle("Critical values y1(y2)")

## -----------------------------------------------------------------------------
## just positive
decision(postmix(treat_prior, n = 24, r = 15), postmix(map, n = 6, r = 3))
## negative
decision(postmix(treat_prior, n = 24, r = 14), postmix(map, n = 6, r = 4))

## -----------------------------------------------------------------------------
r_placebo <- 1
r_treat <- 14

## first obtain posterior distributions...
post_placebo <- postmix(map_robust, r = r_placebo, n = 6)
post_treat <- postmix(treat_prior, r = r_treat, n = 24)

## ...then calculate probability that the difference is smaller than
## zero
prob_smaller <- pmixdiff(post_treat, post_placebo, 0, lower.tail = FALSE)

prob_smaller

prob_smaller > 0.95

## alternativley we can use the decision object
decision(post_treat, post_placebo)

## -----------------------------------------------------------------------------
sessionInfo()

## ----include=FALSE------------------------------------------------------------
options(.user_mc_options)

