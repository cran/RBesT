test_that("write_mix_json and read_mix_json work correctly", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  write_mix_json(nm, tempfile, pretty = TRUE, digits = 1)

  # Check that the file was created
  expect_true(file.exists(tempfile))

  # Deserialize the JSON file back into a mixture object
  mix <- read_mix_json(tempfile)

  # Check that the deserialized object matches the original
  expect_true(all(nm == mix))
})

test_that("write_mix_json and read_mix_json preserve three components", {
  # Create a three-component mixture object
  nm <- mixnorm(
    rob = c(0.2, 0, 2),
    inf = c(0.5, 2, 2),
    weak = c(0.3, -1, 1.5),
    sigma = 5
  )

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  write_mix_json(nm, tempfile, pretty = TRUE, digits = 1)

  # Deserialize the JSON file back into a mixture object
  mix <- read_mix_json(tempfile)

  # Check that the three-component structure survives the round-trip
  expect_identical(ncol(mix), 3L)
  expect_identical(colnames(mix), c("rob", "inf", "weak"))
  expect_true(all(nm == mix))
})

test_that("write_mix_json handles EM objects correctly", {
  # Create a mixture object with EM attributes
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  class(nm) <- c("EM", class(nm)) # Simulate an EM object

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  expect_message(
    write_mix_json(nm, tempfile, pretty = TRUE, digits = 1),
    "Dropping EM information from mixture object before serialization."
  )

  # Deserialize the JSON file back into a mixture object
  mix <- read_mix_json(tempfile)

  # Check that the deserialized object matches the original (without EM attributes)
  expect_true(all(unclass(nm) == mix))
})

test_that("read_mix_json rescaling works correctly", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  nm[1, 1] <- 0.15

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  write_mix_json(nm, tempfile, pretty = TRUE, digits = 2)

  # Deserialize the JSON file with rescaling
  mixrescaled <- read_mix_json(tempfile, rescale = TRUE)

  # Check that the weights sum to 1 after rescaling
  expect_equal(sum(mixrescaled[1, ]), 1)

  # Deserialize the JSON file without rescaling
  mix_no_rescale <- read_mix_json(tempfile, rescale = FALSE)

  # Check that the weights do not sum to 1 without rescaling
  expect_false(sum(mix_no_rescale[1, ]) == 1)
})

test_that("write_mix_json warns whenever a weight is written as zero", {
  ## the robust component weight of 0.03 is rounded to 0 whenever only
  ## a single digit is used for the JSON representation
  bm <- mixbeta(rob = c(0.03, 1, 1), inf = c(0.97, 20, 80))

  tempfile <- withr::local_tempfile(fileext = ".json")

  expect_warning(
    write_mix_json(bm, tempfile, digits = 1),
    "rob"
  )

  ## with sufficient precision no warning must be issued
  expect_no_warning(write_mix_json(bm, tempfile, digits = 5))

  ## weights which jsonlite writes in scientific notation are not
  ## written as zero and must not trigger a warning
  tiny <- mixbeta(rob = c(1E-6, 1, 1), inf = c(1 - 1E-6, 20, 80))
  expect_no_warning(write_mix_json(tiny, tempfile, digits = 4))
  expect_equal(
    unname(read_mix_json(tempfile)[1, ]),
    unname(tiny[1, ]),
    tolerance = 1E-4
  )

  ## negative digits are used by jsonlite as significant digits which
  ## must not be mistaken for a loss of the weights
  expect_no_warning(write_mix_json(bm, tempfile, digits = -1))
})

test_that("read_mix_json loads mixtures with zero weights correctly", {
  bm <- mixbeta(rob = c(0.03, 1, 1), inf = c(0.97, 20, 80))

  tempfile <- withr::local_tempfile(fileext = ".json")

  suppressWarnings(write_mix_json(bm, tempfile, digits = 1))

  expect_silent(mix <- read_mix_json(tempfile))

  ## the mixture must stay intact: no component is dropped, the
  ## weights are finite and sum to unity
  expect_equal(ncol(mix), ncol(bm))
  expect_equal(colnames(mix), colnames(bm))
  expect_true(all(is.finite(mix[1, ])))
  expect_equal(sum(mix[1, ]), 1)
  expect_equal(unname(mix[1, ]), c(0, 1))

  ## and the mixture must be usable
  expect_true(all(is.finite(summary(mix))))
  expect_equal(length(rmix(mix, 5)), 5)
})

test_that("write_mix_json rejects mixtures with all weights rounded to zero", {
  ## 30 equally weighted components each have a weight of 1/30 which
  ## is rounded to 0 whenever only a single digit is used
  Nc <- 30
  bm <- do.call(
    mixbeta,
    lapply(seq_len(Nc), function(i) c(1 / Nc, 10 + i, 80))
  )

  tempfile <- withr::local_tempfile(fileext = ".json")

  expect_error(
    write_mix_json(bm, tempfile, digits = 1),
    "All mixture weights"
  )
})

test_that("read_mix_json rejects mixtures with all weights being zero", {
  bm <- mixbeta(rob = c(0.03, 1, 1), inf = c(0.97, 20, 80))

  tempfile <- withr::local_tempfile(fileext = ".json")

  suppressWarnings(write_mix_json(bm, tempfile, digits = 1))

  ## manually zero out all weights to mimic a mixture written with an
  ## even lower precision
  json <- readLines(tempfile)
  json <- sub('"comp":[[0,1],', '"comp":[[0,0],', json, fixed = TRUE)
  writeLines(json, tempfile)

  expect_error(read_mix_json(tempfile), "All mixture weights")
})

test_that("write_mix_json warns about missing digits argument", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  # Use withr to create a temporary file
  temp_file <- withr::local_tempfile(fileext = ".json")

  # Expect a warning when digits argument is not provided
  expect_warning(
    write_mix_json(nm, temp_file, pretty = TRUE),
    "JSON serialization by default restricts number of digits"
  )
})
