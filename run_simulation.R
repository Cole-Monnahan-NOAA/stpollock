### This files runs the simulation testing component of the analysis.
rm(list=ls())
source("startup.R")
## the simulated data is based on the real data so process that now


## Basic simulation
set.seed(1)
atrend <- c(seq(0,1, len=5), seq(1,-1, len=5))
vtrend <- rep(1,10)
vtrend <- c(4,16,24,16, 18, 18,14, 12,12,12)

## currently depth has no impact but should add that and other covariates later
Nsamples <- 200
st.list <- list(lon=runif(Nsamples,-175, -160), lat=runif(Nsamples, 55,62), beta0=3,
                depth = sample(50:100, size=Nsamples, replace=TRUE))

## ## Run a single replicate in serial for testing
## out <- simulate(replicate=1, st.list=st.list, nyrs=10,
##                 abundance.trend=atrend, plot=TRUE)
## ggplot(out, aes(year, est, color=strata)) + geom_line() + geom_point() +
##   facet_grid(model~space)

chains <- cores <- 6
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
out.parallel <-
  snowfall::sfLapply(1:chains, function(i)
    simulate(replicate=i, st.list=st.list, nyrs=10, abundance.trend=atrend,
             vertical.trend=vtrend, bt.sd=.1, at.sd=.1,
             plot=FALSE))
sfStop()
indices <- do.call(rbind, out.parallel)
indices$group <- with(indices, paste(rep, model, strata, sep='-'))
## indices$year <- indices$year+ ((1:10)/10)[indices$x]
ggplot(indices, aes(year, est, group=group, color=strata)) + geom_line() + geom_jitter() +
  facet_grid(model~space)
