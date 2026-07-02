source_example <- function(
  example,
  env,
  disable_plots = TRUE
) {
  stopifnot(is.environment(env))

  ex_source <- readLines(system.file(
    "examples",
    example,
    package = "RBesT",
    mustWork = TRUE
  ))
  if (disable_plots) {
    ex_source <- grep("plot\\(", ex_source, value = TRUE, invert = TRUE)
  }
  suppressMessages(
    ex <- source(
      textConnection(ex_source),
      local = env,
      echo = FALSE,
      verbose = FALSE
    )
  )
  invisible(ex)
}

run_example_fixture <- function(example, disable_plots = TRUE) {
  env <- new.env(parent = globalenv())
  source_example(
    example = example,
    env = env,
    disable_plots = disable_plots
  )
  as.list(env, all.names = TRUE)
}
