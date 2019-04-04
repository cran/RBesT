#! /bin/bash
#BSUB -J sbc_RBesT
#BSUB -oo sbc_log_out.%J
#BSUB -eo sbc_log_err.%J
#BSUB -q normal
#BSUB -n 1
#BSUB -env "all"
#BSUB -cwd .

export R_HOME=/CHBS/apps/R/3.4.3/lib64/R
export R_ENVIRON_USER=${HOME}/.Renviron_mran
export R_BATCHTOOLS_SEARCH_PATH=${HOME}/batchtools_lsf/

${R_HOME}/bin/Rscript make_reference_rankhist.R
