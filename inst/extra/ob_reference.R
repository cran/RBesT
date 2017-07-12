## Creates OncoBayes reference runs

library(OncoBayes)
library(RBesT)

## source test-gMAP.R but avoid running tests
test_that <<- function(...) {}
context <<- function(...) {}
source("tests/testthat/test-gMAP.R")
rm(context)
rm(test_that)

options(OB.engine="JAGS")

## example 2 from help
OB_normal_data <- data.frame(y=c(0.5,1.5,0.7), s=c(0.3,3.1,0.2), n=c( 10, 34, 51), label=c("s1","s2","s3"))

OB_normal_map <- MAP.Normal(y = OB_normal_data$y,
                            s = OB_normal_data$s,
                            n = OB_normal_data$n,
                            labels = OB_normal_data$label,
                            Dist.tau = "LogNormal",
                            Prior.tau = c(log(0.25), log(2)/1.96),
                            Prior.beta = c(0, 1e2),
                            RndSeed = 10,
                            MCMC = c(2000, 101000, 8, 1),
                            Outfile = "ex2",
                            Plot=FALSE
                            )
sims <- OB_normal_map$R2WinBUGS$sims.list

OB_normal_ref <- lapply(c(sims[c("beta", "tau")], list(theta_resp_pred=sims$theta.pred[,,1])), get_std_quants)

##save(OB_normal_ref, OB_normal_data, file="OncoBayes_ref.rda")

## example 2 from help
OB_binary_map <- MAP.Binary(r = colitis$r,
                            n = colitis$n,
                            study = colitis$study,
                            labels = colitis$study.labels,
                            Dist.tau = "HalfNormal",
                            Prior.tau = c(0, 1),
                            Prior.beta = c(0, 10),
                            RndSeed = 109786,
                            MCMC = c(2000, 101000, 8, 1),
                            Outfile = "ex2",
                            Plot=FALSE
                            )

sims <- OB_binary_map$R2WinBUGS$sims.list

OB_binary_ref <- lapply(c(sims[c("beta", "tau")], list(theta_resp_pred=sims$p.pred[,,1])), get_std_quants)

## example 1 from help
OB_poisson_data <- data.frame(r=c(5,2,7), e=c(10,34,51), label=c("s1", "s2", "s3"))
OB_poisson_map <- MAP.Poisson(r = OB_poisson_data$r,
                              e = OB_poisson_data$e,
                              labels = OB_poisson_data$label,
                              Dist.tau = "LogNormal",
                              Prior.tau = c(log(0.25), log(2)/1.96) ,
                              Prior.beta = c(0, 1e2),
                              RndSeed = 5,
                              REdist =  "normal",
                              Outfile = "ex1",
                              MCMC = c(2000, 101000, 8, 1),
                              Plot=FALSE
                              )
sims <- OB_poisson_map$R2WinBUGS$sims.list
names(sims)

OB_poisson_ref <- lapply(c(sims[c("beta", "tau")], list(theta_resp_pred=sims$lambda.pred[,,1])), get_std_quants)

## only keep the reference objects and data
rm(sims, OB_poisson_map, OB_binary_map, OB_normal_map)

OB_session <- sessionInfo()
