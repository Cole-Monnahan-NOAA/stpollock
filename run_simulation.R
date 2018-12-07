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
atrend <- c(seq(0,1, len=5), seq(1,-1, len=5))
## currently depth has no impact but should add that and other covariates later
Nsamples <- 200
st.list <- list(lon=runif(Nsamples,-175, -160), lat=runif(Nsamples, 55,62), beta0=3,
                depth = sample(50:100, size=Nsamples, replace=TRUE))
out <- simulate(replicate=1, st.list=st.list, nyrs=10,
                abundance.trend=atrend, plot=FALSE)

par(mfrow=c(1,2))
plot(1:10, atrend)
plot(1:10, out$vast.full$index$value, ylim=c(0, 6e6))
with(out$vast.full$index, lines(1:10, value+1.96*se))
with(out$vast.full$index, plot(1:10, value-1.96*se))

chains <- cores <- 6
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
mcmc.out <- snowfall::sfLapply(1:chains, function(i)
 simulate(replicate=i, st.list=st.list, nyrs=10, abundance.trend=atrend))
sfStop()

