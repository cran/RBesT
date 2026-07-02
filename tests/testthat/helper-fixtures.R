#' Resolve A Cached MCMC gMAP Test Fixture Path
#'
#' Fixture sources live in `tests/testthat/fixtures-mcmc-src/` and are built by
#' `make test-fixtures` into `tests/testthat/fixtures-mcmc/`. Generated `.rds`
#' files are local cache artifacts and are not committed.
#'
#' @param name Fixture name without the `.rds` suffix.
#' @return Absolute path to the cached fixture file.
mcmc_gmap_fixture_path <- function(name) {
  testthat::test_path("fixtures-mcmc", paste0(name, ".rds"))
}

#' Load A gMAP Test Fixture
#'
#' Tests must explicitly request whether they need a compact fixture or a real
#' cached MCMC fixture. Compact fixtures are intended to be cheap and CRAN-safe;
#' MCMC fixtures remain non-CRAN local cache artifacts.
#'
#' @param name Fixture name without the generated suffix.
#' @param type Fixture backend, either `"compact"` or `"mcmc"`.
#' @return A validated `gMAP` object with stored posterior draws.
load_gmap_fixture <- function(name, type) {
  type <- match.arg(type, c("compact", "mcmc"))

  switch(
    type,
    compact = load_compact_gmap_fixture(name),
    mcmc = load_mcmc_gmap_fixture(name)
  )
}

#' Load A Cached MCMC gMAP Test Fixture
#'
#' Fixture-backed tests are treated as non-CRAN tests because the generated
#' fixture cache is intentionally not shipped with the package source, saving
#' source-package storage size. If the local cache is missing, the test is
#' skipped with instructions to run `make -j4 test-fixtures`.
#'
#' Cached fixture `.rds` files are trusted local build artifacts generated from
#' committed recipes under `tests/testthat/fixtures-mcmc-src/`. They are ignored by
#' git and excluded from source builds, so fixture-backed tests should not load
#' `.rds` files obtained from untrusted sources.
#'
#' @param name Fixture name without the `.rds` suffix.
#' @return A validated `gMAP` object with stored posterior draws.
load_mcmc_gmap_fixture <- function(name) {
  testthat::skip_on_cran()

  path <- mcmc_gmap_fixture_path(name)
  if (!file.exists(path)) {
    testthat::skip(
      paste0(
        "Missing gMAP fixture: ", path, "\n",
        "Run `make -j4 test-fixtures` first to create cached test fixtures."
      )
    )
  }

  fixture <- readRDS(
    path,
    refhook = function(...) {
      stop("External references are not allowed in test fixtures.", call. = FALSE)
    }
  )
  if (!inherits(fixture, "gMAP")) {
    stop("Fixture is not a gMAP object: ", path, call. = FALSE)
  }
  if (is.null(fixture$draws)) {
    stop("Fixture does not contain posterior draws: ", path, call. = FALSE)
  }

  fixture
}
