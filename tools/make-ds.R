#! /usr/bin/env Rscript

rm(list=ls())

## generate example data sets
make_ds <- function() {
    colitis <<- data.frame(study=c("Van_assche", "Feagan", "Rutgeerts-1", "Rutgeerts-2"),
                           n = c(56,63,121,123),
                           r = c(6,9,18,7)
                           )

    AS <<- data.frame(study=paste("Study", 1:8),
                      n=c(107,44,51,39,139,20,78,35),
                      r=c(23,12,19,9,39,6,9,10))
    
    transplant <<- data.frame(study=paste("Study", 1:11),
                              n=c( 33, 45, 74,103,140, 49, 83, 59, 22,109,213),
                              r=c(  6,  8, 17, 28, 26,  8, 22,  8,  6, 16, 53))

    crohn <<- dat <- data.frame(study=c("Gastr06","AIMed07","NEJM07","Gastr01a","APhTh04","Gastr01b"),
                                n=c(74, 166, 328, 20, 25, 58),
                                y=c(-51, -49, -36, -47, -90, -54))

    ## use_data expects it's data sets in the global env (which is why
    ## we do <<-)
    use_data(AS, transplant, colitis, crohn, overwrite=TRUE)
}


make_internal_ds <- function() {
    if(!file.exists("inst/extra/ob_reference.rda"))
        make_OB_reference()
    load("inst/extra/ob_reference.rda",envir=.GlobalEnv)

    use_data(
        OB_binary_ref,
        OB_normal_data,
        OB_normal_ref,
        OB_poisson_data,
        OB_poisson_ref,
        OB_session,
        internal=TRUE, overwrite=TRUE)
}

## create reference OncoBayes runs
make_OB_reference <- function() {
    ob <- new.env()
    source("inst/extra/ob_reference.R", local=ob)
    save(list=ls(ob), envir=ob, file="inst/extra/ob_reference.rda")
}

library(devtools)

## cleanup first
if(file.exists("R/sysdata.rda"))
    file.remove("R/sysdata.rda")

for(rda in dir("data", pattern="*rda", full.names=TRUE))
    file.remove(rda)

##load_all()

make_ds()
##make_OB_reference()
make_internal_ds()
