## test the analytical OC function via brute force simulation

## expect results to be 1% exact
pos1S_tolerance <- function() 1e-2

## we test here that the PoS indeed averages over the predictive of an
## informative prior the conditional power.

pos1S_sample_size <- function() 1E4

pos1S_reference_scenario <- function() {
  s <- 2
  theta_ni <- 0.4
  alpha <- 0.05
  n1 <- 155
  prior <- mixnorm(c(1, 0, 100), sigma = s)

  list(
    decA = decision1S(1 - alpha, theta_ni, lower.tail = TRUE),
    decAU = decision1S(1 - alpha, theta_ni, lower.tail = FALSE),
    ia_dist = postmix(prior, m = log(0.9), se = s / sqrt(40)),
    n1 = n1,
    prior = prior,
    s = s
  )
}

pos1S_beta_scenario <- function() {
  sc <- pos1S_reference_scenario()
  prior <- mixbeta(c(1, 1, 1))
  list(
    decA = sc$decA,
    decAU = sc$decAU,
    ia_dist = postmix(prior, r = 20, n = 50),
    n1 = sc$n1,
    prior = prior
  )
}

pos1S_gamma_scenario <- function() {
  alpha <- 0.05
  prior <- mixgamma(c(1, 1, 1), param = "mn")
  list(
    dec_count = decision1S(1 - alpha, 1, lower.tail = TRUE),
    dec_countU = decision1S(1 - alpha, 1, lower.tail = FALSE),
    ia_dist = postmix(prior, m = 0.9, n = 40),
    n1 = 155,
    prior = prior
  )
}

pos1S_extreme_normal_scenario <- function() {
  sc <- pos1S_reference_scenario()
  prior <- mixnorm(c(1, 0, 1), sigma = sc$s, param = "mn")
  list(
    decA = sc$decA,
    decAU = sc$decAU,
    ia_dist = postmix(prior, m = -1, n = 162),
    n = 459 - 162,
    prior = prior
  )
}

test_pos1S <- function(prior, ia_dist, n, dec, decU) {
  withr::local_seed(873621)

  ## the PoS is the expected value of the condition power integrated
  ## over the interim density which is what we check here
  cpo_analytic <- oc1S(prior, n, dec)
  pos_analytic <- pos1S(prior, n, dec)
  samp <- rmix(ia_dist, pos1S_sample_size())
  pos_mc <- mean(cpo_analytic(samp))
  expect_true(all(abs(pos_mc - pos_analytic(ia_dist)) < pos1S_tolerance()))
  lower.tail <- attr(dec, "lower.tail")
  if (lower.tail) {
    test_pos1S(prior, ia_dist, n, decU)
  }
}


test_that("Normal PoS 1 sample function matches MC integration of CPO", {
  sc <- pos1S_reference_scenario()
  test_pos1S(sc$prior, sc$ia_dist, sc$n1, sc$decA, sc$decAU)
})

test_that("Binomial PoS 1 sample function matches MC integration of CPO", {
  sc <- pos1S_beta_scenario()
  test_pos1S(sc$prior, sc$ia_dist, sc$n1, sc$decA, sc$decAU)
})

test_that("Poisson PoS 1 sample function matches MC integration of CPO", {
  sc <- pos1S_gamma_scenario()
  test_pos1S(sc$prior, sc$ia_dist, sc$n1, sc$dec_count, sc$dec_countU)
})

test_that("Normal PoS 1 sample function matches MC integration of CPO (more extreme case)", {
  sc <- pos1S_extreme_normal_scenario()
  test_pos1S(sc$prior, sc$ia_dist, sc$n, sc$decA, sc$decAU)
})

test_that("Mixed lower.tail usage works for normal PoS calculation", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  post_ia <- postmix(prior, m = -1, n = 15)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for binomial PoS calculation", {
  prior <- mixbeta(rob = c(0.2, 0.5, 0.5), inf = c(0.8, 0.5, 0.5))
  post_ia <- postmix(prior, r = 20, n = 50)

  dec_lower <- decision1S(pc = 0.5, qc = 0.8, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(0.8, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for Poisson PoS calculation", {
  prior <- mixgamma(
    rob = c(0.2, 0.5, 0.5),
    inf = c(0.8, 0.5, 0.5),
    param = "mn"
  )
  post_ia <- postmix(prior, m = 0.9, n = 40)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})
