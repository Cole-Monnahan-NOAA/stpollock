rm(list=ls())
library(TMB)
library(VAST)
library(reshape2)
library(ggplot2)
library(TMBhelper)
library(snowfall)
source("simulator.R")

## Basic simulation
set.seed(1)
atrend <- sort(rnorm(10))
st.list <- list(lon=runif(100), lat=runif(100), beta0=3,
                depth = sample(50:100, size=100, replace=TRUE))

out <- simulate(replicate=1, st.list=st.list, nyrs=10, abundance.trend=atrend)
plot(atrend,out[[1]]$index)


chains <- cores <- 3
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
mcmc.out <- snowfall::sfLapply(1:chains, function(i)
 simulate(replicate=i, st.list=st.list, nyrs=10, abundance.trend=atrend))
sfStop()
