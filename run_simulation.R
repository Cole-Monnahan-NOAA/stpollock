rm(list=ls())
library(TMB)
library(VAST)
library(reshape2)
library(ggplot2)
library(TMBhelper)
library(snowfall)
library(maps)
library(mapdata)
Version <- "VAST_v5_3_0"
source("simulator.R")

## Basic simulation
set.seed(1)
atrend <- rnorm(10)
## currently depth has no impact but should add that and other covariates later
st.list <- list(lon=runif(100,-175, -160), lat=runif(100, 55,62), beta0=3,
                depth = sample(50:100, size=100, replace=TRUE))
out <- simulate(replicate=1, st.list=st.list, nyrs=10, abundance.trend=atrend)
plot(atrend,out[[1]]$index)
tmp-out$vast.full$Opt$par

chains <- cores <- 6
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
mcmc.out <- snowfall::sfLapply(1:chains, function(i)
 simulate(replicate=i, st.list=st.list, nyrs=10, abundance.trend=atrend))
sfStop()

