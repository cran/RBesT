library(pkgbuild)
## ensure that the current dev version of RBesT is loaded
pkgbuild::compile_dll("../..")

library(assertthat)
library(batchtools)
source("sbc_tools.R")
set.seed(453453)

#' according to the docs this speeds up the reduce step
options(batchtools.progress = FALSE)

#' Evaluate dense and sparse data-scenario

#' Dense: 10 trials with 40 entries each
dense_data  <- list(group=rep(1:10, each=40))
#' Sparse: 2 trials with 40 entries each
sparse_data <- list(group=rep(1:2,  each=40))

reg <- makeExperimentRegistry(file.dir = tempfile("sbc_"),
                              ## use the default configured batchtools configuration batchtools.conf
                              ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
                              ## - current working directory
                              ## - $(HOME)/.batchtools.conf
                              seed = 47845854,
                              ## our worker functions and package loading
                              source="sbc_tools.R"
                              )

if(FALSE) {
    ## for debugging here
    removeProblems("dense")
    removeProblems("sparse")
}

addProblem("dense",
           data=dense_data,
           fun=simulate_fake,
           seed=2345,
           ## caching speeds up the reduce step
           cache=TRUE
           )

addProblem("sparse",
           data=sparse_data,
           fun=simulate_fake,
           seed=2345,
           ## caching speeds up the reduce step
           cache=TRUE
           )

addAlgorithm("RBesT", fit_rbest)

## family, mean_mu, sd_mu, sd_tau, samp_sd
scenarios <- data.frame(
    family=c("binomial", "gaussian", "poisson"),
    mean_mu=c(-1, 0, 0),
    sd_mu=c(1),
    sd_tau=c(rep(0.5, 3), rep(1, 3)),
    samp_sd=c(1),
    stringsAsFactors=FALSE)

pdes <- list(sparse=scenarios, dense=scenarios)
ades <- list(RBesT=data.frame())

#' Add the defined problems and analysis methods to the registry and
#' set the number of replications:
S  <- 1E4L
addExperiments(pdes, ades, repls=S)

summarizeExperiments()

if(FALSE) {
    ## used for debugging
    job1 <- testJob(1)
    job2 <- testJob(6)
    job3 <- testJob(11)

    job <- makeJob(1)

    fit_rbest(problem_data_dense, job, job$instance )
}


#'
#' Chunk the jobs into 600 chunks to run
#'
ids <- getJobTable()
ids[, chunk:=chunk(job.id, 600)]

#' Once things run fine let's submit this work to the cluster.
submitJobs(ids)

#' Wait for results.
waitForJobs()

#' Check status:
getStatus()

#' Ensure that no error occured
assert_that(nrow(findErrors()) == 0)

#' Collect results.
calibration_data <- ijoin(
    ## grab job parameters
    unwrap(getJobPars()),
    unwrap(reduceResultsDataTable(fun=function(x) c(rank=as.list(x$rank),
                                                    list(n_divergent=x$n_divergent,
                                                         min_Neff=ceiling(x$min_Neff))) ))
)

calibration_data[,algorithm:=NULL]

#' Bin raw data as used in the analysis.
scale64  <- scale_ranks(1024, 2^4)
B <- 1024L / 2^4
calibration_data_binned <- calibration_data[, scale64(.SD), by=c("problem", "family", "sd_tau")]

#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_data_binned <- as.data.frame(calibration_data_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout=TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz="UTC")

calibration <- list(raw=calibration_data,
                    data=calibration_data_binned,
                    S=S,
                    B=B,
                    git_hash=git_hash,
                    created=created)

saveRDS(calibration, file="calibration.rds")

library(tools)
md5 <- md5sum("calibration.rds")
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"), file="calibration.md5")

#' Cleanup
removeRegistry(0)

#' Session info
sessionInfo()
