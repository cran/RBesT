#!/bin/bash -x

[ -z "$R_HOME" ] && R_HOME=`R RHOME`

## ensure that NAMESPACE contains load directive
echo -e "# Generated by roxygen2: do not edit by hand\nuseDynLib(RBesT, .registration = TRUE)\n" > NAMESPACE

## compile RBesT dll
"${R_HOME}/bin/R" --slave -e 'library(pkgbuild); pkgbuild::compile_dll()'

## create internal data-sets
"${R_HOME}/bin/R" --slave --file=tools/make-ds.R

## create documentation Rd files
"${R_HOME}/bin/R" --slave -e 'library(roxygen2); roxygen2::roxygenize()'

## create SBC report
"${R_HOME}/bin/R" --slave -e 'library(rmarkdown); setwd("inst/sbc/"); rmarkdown::render("sbc_report.R")'

