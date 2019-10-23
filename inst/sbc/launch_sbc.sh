#! /bin/bash
#BSUB -J sbc_RBesT
#BSUB -oo sbc_log_out.%J
#BSUB -eo sbc_log_err.%J
#BSUB -n 1
## execution time of up to 8h (short queue) and 5G RAM
#BSUB -q short
#BSUB -M 5000
#BSUB -env "all"
#BSUB -cwd .

export R_HOME=/CHBS/apps/R/3.4.3/lib64/R
export R_ENVIRON_USER=${HOME}/.Renviron_mran
export R_BATCHTOOLS_SEARCH_PATH=${HOME}/batchtools_lsf/

## todo: script which starts the jobs with few ressources and then a
## script which collects results => then both of these can run in the
## short queue and we do not need to hold ressources which simply wait

${R_HOME}/bin/Rscript make_reference_rankhist.R
