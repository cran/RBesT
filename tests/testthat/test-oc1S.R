## test the analytical OC function via brute force simulation

## expect results to be 1% exact
oc1S_tolerance <- function() 1e-2

## here we test against the reference Neuenschwander et al.,
## Statist. Med. 2011, 30, 1618 (*the* double criterion paper)

oc1S_reference_scenario <- function() {
  s <- 2
  theta_ni <- 0.4
  theta_a <- 0
  alpha <- 0.05
  beta <- 0.2
  n1 <- 155
  c1 <- theta_ni - qnorm(1 - alpha) * s / sqrt(n1)

  ## standard NI design, tests only statistical significance to be
  ## smaller than theta_ni with 1-alpha certainty
  decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
  prior <- mixnorm(c(1, 0, 100), sigma = s)

  ## double criterion with indecision point equal to the classical boundary
  theta_c <- c1
  dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
  dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
  decComb <- decision1S(
    c(1 - alpha, 0.5),
    c(theta_ni, theta_c),
    lower.tail = TRUE
  )

  list(
    alpha = alpha,
    beta = beta,
    c1 = c1,
    dec1 = dec1,
    dec1b = decision1S(1 - alpha, theta_ni, lower.tail = FALSE),
    dec2 = dec2,
    decA = decA,
    decComb = decComb,
    n1 = n1,
    prior = prior,
    thetaA = c(theta_a, theta_ni),
    thetaD = c(theta_c, theta_ni),
    theta_ni = theta_ni
  )
}

test_scenario <- function(oc_res, ref) {
  resA <- oc_res - ref
  expect_true(all(abs(resA) < oc1S_tolerance()))
}

test_that("Classical NI design critical value", {
  sc <- oc1S_reference_scenario()
  expect_true(
    abs(decision1S_boundary(sc$prior, sc$n1, sc$decA) - sc$c1) <
      oc1S_tolerance()
  )
})

## n set to give power 80% to detect 0 and type I error 5% for no
## better than theta_ni
test_that("Classical NI design at target    sample size for OCs", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, sc$n1, sc$decA)(sc$thetaA), c(1 - sc$beta, sc$alpha))
})
test_that("Classical NI design at increased sample size for OCs", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 233, sc$decA)(sc$thetaA), c(1 - 0.08, sc$alpha))
})
test_that("Classical NI design at decreased sample size for OCs", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 77, sc$decA)(sc$thetaA), c(1 - 0.45, sc$alpha))
})

## since theta_c == c1, both decision criteria are the same for n =
## 155
test_that("Double criterion NI design at target    sample size for OCs, combined      ", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, sc$n1, sc$decComb)(sc$thetaD), c(0.50, sc$alpha))
})
test_that("Double criterion NI design at target    sample size for OCs, stat criterion", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, sc$n1, sc$dec1)(sc$thetaD), c(0.50, sc$alpha))
})
test_that("Double criterion NI design at target    sample size for OCs, mean criterion", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, sc$n1, sc$dec2)(sc$thetaD), c(0.50, sc$alpha))
})

## at an increased sample size only the mean criterion is active
test_that("Double criterion NI design at increased sample size for OCs, combined      ", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 233, sc$decComb)(sc$thetaD), c(0.50, 0.02))
})
test_that("Double criterion NI design at increased sample size for OCs, mean criterion", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 233, sc$dec2)(sc$thetaD), c(0.50, 0.02))
})

## at a  decreased sample size only the stat criterion is active
test_that("Double criterion NI design at decreased sample size for OCs, combined      ", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 78, sc$decComb)(sc$thetaD), c(1 - 0.68, sc$alpha))
})
test_that("Double criterion NI design at decreased sample size for OCs, stat criterion", {
  sc <- oc1S_reference_scenario()
  test_scenario(oc1S(sc$prior, 78, sc$dec1)(sc$thetaD), c(1 - 0.68, sc$alpha))
})

## design object, decision function, posterior function must return
## posterior after updatding the prior with the given value
test_critical_discrete <- function(crit, decision, posterior) {
  lower.tail <- attr(decision, "lower.tail")
  if (lower.tail) {
    expect_equal(decision(posterior(crit - 1)), 1)
    expect_equal(decision(posterior(crit)), 1)
    expect_equal(decision(posterior(crit + 1)), 0)
  } else {
    expect_equal(decision(posterior(crit - 1)), 0)
    expect_equal(decision(posterior(crit)), 0)
    expect_equal(decision(posterior(crit + 1)), 1)
  }
}

oc1S_binary_scenario <- function() {
  sc <- oc1S_reference_scenario()
  beta_prior <- mixbeta(c(1, 1, 1))
  list(
    alpha = sc$alpha,
    crit_lower = decision1S_boundary(beta_prior, 1000, sc$dec1),
    crit_upper = decision1S_boundary(beta_prior, 1000, sc$dec1b),
    dec_lower = sc$dec1,
    dec_upper = sc$dec1b,
    design_lower = oc1S(beta_prior, 1000, sc$dec1),
    design_upper = oc1S(beta_prior, 1000, sc$dec1b),
    posterior = function(r) postmix(beta_prior, r = r, n = 1000),
    theta_ni = sc$theta_ni
  )
}

test_that("Binary type I error rate", {
  sc <- oc1S_binary_scenario()
  test_scenario(sc$design_lower(sc$theta_ni), sc$alpha)
})
test_that("Binary crticial value, lower.tail=TRUE", {
  sc <- oc1S_binary_scenario()
  test_critical_discrete(sc$crit_lower, sc$dec_lower, sc$posterior)
})
test_that("Binary crticial value, lower.tail=FALSE", {
  sc <- oc1S_binary_scenario()
  test_critical_discrete(sc$crit_upper, sc$dec_upper, sc$posterior)
})

test_that("Binary boundary case, lower.tail=TRUE", {
  sc <- oc1S_binary_scenario()
  expect_numeric(
    sc$design_lower(1),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})
test_that("Binary boundary case, lower.tail=FALSE", {
  sc <- oc1S_binary_scenario()
  expect_numeric(
    sc$design_upper(0),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})

oc1S_poisson_scenario <- function() {
  sc <- oc1S_reference_scenario()
  gamma_prior <- mixgamma(c(1, 2, 2))
  dec_count <- decision1S(1 - sc$alpha, 1, lower.tail = TRUE)
  dec_countB <- decision1S(1 - sc$alpha, 1, lower.tail = FALSE)
  list(
    alpha = sc$alpha,
    crit_lower = decision1S_boundary(gamma_prior, 1000, dec_count),
    crit_upper = decision1S_boundary(gamma_prior, 1000, dec_countB),
    dec_lower = dec_count,
    dec_upper = dec_countB,
    design_lower = oc1S(gamma_prior, 1000, dec_count),
    posterior = function(m) postmix(gamma_prior, m = m / 1000, n = 1000)
  )
}

test_that("Poisson type I error rate", {
  sc <- oc1S_poisson_scenario()
  test_scenario(sc$design_lower(1), sc$alpha)
})
test_that("Poisson critical value, lower.tail=TRUE", {
  sc <- oc1S_poisson_scenario()
  test_critical_discrete(sc$crit_lower, sc$dec_lower, sc$posterior)
})
test_that("Poisson critical value, lower.tail=FALSE", {
  sc <- oc1S_poisson_scenario()
  test_critical_discrete(sc$crit_upper, sc$dec_upper, sc$posterior)
})

test_that("Mixed lower.tail usage works for normal OC calculation", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  result_lower <- oc1S(
    prior,
    n = 50,
    decision = dec_lower
  )

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  result_upper <- oc1S(
    prior,
    n = 50,
    decision = dec_upper
  )

  decMixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  result <- oc1S(prior, 50, decMixed)

  theta_grid <- seq(-5, 5, length.out = 50)
  vals_lower <- result_lower(theta_grid)
  vals_upper <- result_upper(theta_grid)
  vals_mixed <- result(theta_grid)

  expected_mixed <- vals_lower - (1 - vals_upper)
  expect_equal(vals_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for binomial OC calculation", {
  prior <- mixbeta(rob = c(0.2, 2, 2), inf = c(0.8, 5, 5))

  dec_lower <- decision1S(pc = 0.5, qc = 0.7, lower.tail = TRUE)
  result_lower <- oc1S(
    prior,
    n = 50,
    decision = dec_lower
  )

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  result_upper <- oc1S(
    prior,
    n = 50,
    decision = dec_upper
  )

  decMixed <- decision1S(
    qc = c(0.7, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  result <- oc1S(prior, 50, decMixed)

  theta_grid <- seq(0, 1, length.out = 50)
  vals_lower <- result_lower(theta_grid)
  vals_upper <- result_upper(theta_grid)
  vals_mixed <- result(theta_grid)

  expected_mixed <- vals_lower - (1 - vals_upper)
  expect_equal(vals_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for Poisson OC calculation", {
  prior <- mixgamma(rob = c(0.2, 2, 2), inf = c(0.8, 5, 5))

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  result_lower <- oc1S(
    prior,
    n = 50,
    decision = dec_lower
  )

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  result_upper <- oc1S(
    prior,
    n = 50,
    decision = dec_upper
  )

  decMixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  result <- oc1S(prior, 50, decMixed)

  theta_grid <- seq(0, 10, length.out = 50)
  vals_lower <- result_lower(theta_grid)
  vals_upper <- result_upper(theta_grid)
  vals_mixed <- result(theta_grid)

  expected_mixed <- vals_lower - (1 - vals_upper)
  expect_equal(vals_mixed, expected_mixed)
})
