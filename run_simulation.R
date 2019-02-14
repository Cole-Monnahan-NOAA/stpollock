### This files runs the simulation testing component of the analysis.
rm(list=ls())
source("startup.R")
## the simulated data is based on the real data so process that now


## Basic simulation
set.seed(1)
atrend <- c(seq(0,1, len=5), seq(1,-1, len=5))
vtrend <- c(4,16,24,16, 20, 25,14, 12,12,12)
vtrend <- c(rep(5, 5), rep(5,5)); atrend <- rep(1, 10)

## currently depth has no impact but should add that and other covariates later
Nsamples <- 1000
st.list <-
  list(lon=runif(Nsamples,-175, -160), lat=runif(Nsamples, 55,62),
       beta0=5, depth = sample(50:100, size=Nsamples, replace=TRUE),
       sd.process=.005, nyrs=10,  bt.sd=.001, at.sd=.001, n_x=10)

## ## Run a single replicate in serial for testing
out <- simulate(replicate=12, st.list=st.list,
                atrend=atrend, vtrend=vtrend, plot=TRUE)
## ggplot(out, aes(year, est, color=strata)) + geom_line() + geom_point() +
##   facet_grid(model~space)

cores <- 6
chains <- cores
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
out.parallel <-
  ##sapply(1:chains, function(i)
  snowfall::sfLapply(1:chains, function(i)
    simulate(replicate=i, st.list=st.list, atrend=atrend,
             vtrend=vtrend, plot=TRUE))
sfStop()
indices <- do.call(rbind, out.parallel)

## indices$year <- indices$year+ ((1:10)/10)[indices$x]
ggplot(indices, aes(year, est, group=group, color=model)) + geom_line()+
  facet_grid(strata~space)

## have to manipulate this one a bit
x <- droplevels(subset(indices, model == 'combined'))
y <- melt(x, id.vars=c('group2', 'year', 'strata2'), measure.vars=c('strata.est', 'truth'))
y$group3 <- paste0(y$group2, '_', y$variable)
ggplot(y, aes(year, value, group=group3, color=variable)) + geom_line() +
  facet_wrap('strata2')
