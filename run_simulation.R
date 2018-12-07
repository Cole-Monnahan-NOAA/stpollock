rm(list=ls())
library(TMB)
library(VAST)
library(reshape2)
library(ggplot2)
library(plyr)
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
                abundance.trend=atrend, plot=TRUE)

ggplot(out$vast.full$index, aes(year, y=value)) +
  geom_errorbar(aes(ymin=value-1.96*se, ymax=value+1.96*se))

chains <- cores <- 6
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
out.parallel <- snowfall::sfLapply(1:chains, function(i)
 simulate(replicate=i, st.list=st.list, nyrs=10, abundance.trend=atrend))
sfStop()
indices <- ldply(1:cores, function(x) cbind(x,out.parallel[[x]]$vast.full$index))
indices$year <- indices$year+ ((1:10)/10)[indices$x]
ggplot(indices, aes(year, y=value, group=x)) +
  geom_errorbar(aes(ymin=value-1.96*se, ymax=value+1.96*se))
