### This files runs the simulation testing component of the analysis.

rm(list=ls())
source("startup.R")
chains <- 4
options(mc.cores = chains)


## We condition the OM on the fitted base case model
load('fit_basecase_100/Save.RData')
load('fit_basecase_100/Record.RData')
## Rebuild the Obj with the MLEs
mle <- Save$est$est
par.names <- as.character(Save$est$par)
control <- Record$control; control$make_plots <- FALSE
savedir <- paste0(getwd(), '/simulations/simfit_OM')
source("prepare_inputs.R")
Obj$par <- mle
Obj$fn(mle) -  Save$Opt$objective ## check MLE fit
## Obj$gr(mle) ## check that the gradients are 0
pars.all <- Obj$env$last.par ## the full parameter vector

## Build a simulated OM by changing the betas and zero-centering the
## spatial and spatiotemporal effects
Omegainput1_sf <- (pars.all[grep('Omegainput1_sf', names(pars.all))])
tmp <- matrix(Omegainput1_sf, ncol=3, byrow=FALSE)
tmp <- apply(tmp, 2, function(x) x-mean(x))
Omegainput1_sf_centered <- as.numeric(tmp)
pars.all[grep('Omegainput1_sf', names(pars.all))] <-
  Omegainput1_sf_centered
Omegainput2_sf <- (pars.all[grep('Omegainput2_sf', names(pars.all))])
tmp <- matrix(Omegainput2_sf, ncol=3, byrow=FALSE)
tmp <- apply(tmp, 2, function(x) x-mean(x))
Omegainput2_sf_centered <- as.numeric(tmp)
pars.all[grep('Omegainput2_sf', names(pars.all))] <-
  Omegainput2_sf_centered
Epsiloninput1_sft <- (pars.all[grep('Epsiloninput1_sft', names(pars.all))])
tmp <- array(Epsiloninput1_sft, dim=c(116, 3, 12))
tmp <- apply(tmp, c(2,3), function(x) x-mean(x))
Epsiloninput1_sft_centered <- as.vector(tmp)
pars.all[grep('Epsiloninput1_sft', names(pars.all))] <-
  Epsiloninput1_sft_centered
Epsiloninput2_sft <- (pars.all[grep('Epsiloninput2_sft', names(pars.all))])
tmp <- array(Epsiloninput2_sft, dim=c(116, 3, 12))
tmp <- apply(tmp, c(2,3), function(x) x-mean(x))
Epsiloninput2_sft_centered <- as.vector(tmp)
pars.all[grep('Epsiloninput2_sft', names(pars.all))] <-
  Epsiloninput2_sft_centered
## now the random effects should be centered so the betas drive the trends
## completely (right?)

beta1_ft <- (pars.all[grep('beta1_ft', names(pars.all))])
beta1_ft <- matrix(beta1_ft, ncol=3, byrow=TRUE)
beta1 <- cbind(c(seq(1,-2.5, len=12)),
               c(seq(0, 3.5, len=12)),
               c(seq(.5,3, len=12)))
pars.all[grep('beta1_ft', names(pars.all))] <- as.vector(t(beta1))
## beta2_ft <- (pars.all[grep('beta2_ft', names(pars.all))])
## beta2_ft <- matrix(beta2_ft, ncol=3, byrow=TRUE)
## beta2 <- cbind(seq(0,.5, len=12), seq(.5, 3, len=12), seq(0, 1, len=12))
## pars.all[grep('beta2_ft', names(pars.all))] <- as.vector(t(beta2))
## par(mfrow=c(2,2))
## matplot(beta1)
## matplot(beta1_ft)
## matplot(beta2)
## matplot(beta2_ft)

## ## Put data into the ATS for missing years so it's sampled in the
## ## simulation !!! this actually breaks it b/c the mesh construction is
## ## different so the random effects are in the wrong order
## source('load_data.R')
## tmp2 <- subset(Data_Geostat, Gear == 'Acoustic_3-16' & Year %in% c(2010, 2012, 2014, 2016))
## tmp2$Year <- tmp2$Year+1
## DF2 <- rbind(DF2, tmp2)
## tmp3 <- subset(Data_Geostat, Gear == 'Acoustic_16-surface' & Year %in% c(2010, 2012, 2014, 2016))
## tmp3$Year <- tmp3$Year+1
## DF3 <- rbind(DF3, tmp3)
## ## Have to rebuild the object since the data changed
## control$simdata <- TRUE
## source("prepare_inputs.R")
## table(Data_Geostat$Gear) ## check the simdata was used
sim <- Obj$report(pars.all)
index.sim <- t(sim$Index_cyl[,,1])
par(mfrow=c(1,2))
matplot(log(cbind(index.sim, rowSums(index.sim))), ylab='log index')
matplot(t(apply(index.sim, 1, function(x) cumsum(x/sum(x)))), type='l', ylim=c(0,1))
## Now when I call a new vector it'll have expected values in the missing
## years and thus be sampled from them.

Data_Geostat0 <- Data_Geostat
### Loop through replicates of the simulation
for(i in 2:5){
## Simulate the sampling process.
set.seed(i)
encounters <- rbinom(length(sim$R1_i), size=1, prob=sim$R1_i)
## Carefully index to get the right SigmaObs by gear type
sigma.obs <- sim$SigmaM[,1][as.numeric(Data_Geostat0$Gear)]
logcatches <- rnorm(length(sim$R2_i), log(sim$R2_i)-sigma.obs^2/2, sigma.obs)
catches <- exp(logcatches)*encounters

## Rebuild the Obj with the new simulated data
Data_Geostat <- Data_Geostat0
Data_Geostat$Catch_KG <- catches
DF1 <- subset(Data_Geostat, Gear=='Trawl')
DF2 <- subset(Data_Geostat, Gear=='Acoustic_3-16')
DF3 <- subset(Data_Geostat, Gear=='Acoustic_16-surface')
DF1$knot_i <- DF2$knot_i <- DF3$knot_i <- NULL

## Master settings which match the fitted OM
getsd <- FALSE
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='combined',
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
control$simdata <- TRUE
savedir <- paste0(getwd(), '/simulations/simfit_', i, '_combined')
source("prepare_inputs.R")
saveRDS(sim, paste0(savedir, '/simreport.RDS'))
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Repeat with just the BTS
control$model <- 'bts'; control$make_plots <- FALSE
savedir <- paste0(getwd(), '/simulations/simfit_', i, '_bts')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)
saveRDS(sim, paste0(savedir, '/simreport.RDS'))

## Repeat with just the ATS
control$model <- 'ats'
savedir <- paste0(getwd(), '/simulations/simfit_', i, '_ats')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)
saveRDS(sim, paste0(savedir, '/simreport.RDS'))
} ## end of loop over the replicates


## Process and compare the fits to the simulated truth
index.sim <- sim$ln_ColeIndex_cy
index <- t(sim$Index_cyl[,,1])
matplot(log(index))
## load('simulations/simfit_1_combined/Save.RData')
## indexc <- t(Save$Report$Index_cyl[,,1])
## load('simulations/simfit_1_bts/Save.RData')
## indexb <- Save$Report$Index_cyl[,,1]
## load('simulations/simfit_1_ats/Save.RData')
## indexa <- Save$Report$Index_cyl[,,1]

library(magrittr)
library(tidyr)
library(dplyr)
indexc <- readRDS('simulations/simfit_1_combined/index.mcmc.RDS')
indexb <- readRDS('simulations/simfit_1_bts/index.mcmc.RDS')
indexa <- readRDS('simulations/simfit_1_ats/index.mcmc.RDS')
bts.sim <- log(rowSums(index[, 1:2]))
bts.est <- indexb$index.strata$est
bts.est2 <- indexc$index.gear %>% filter(gear=='BTS') %>% select(est) %>% .$est
ats.sim <- log(rowSums(index[, 2:3]))
ats.est <- indexa$index.strata$est
ats.est2 <- indexc$index.gear %>% filter(gear=='ATS') %>% select(est) %>% .$est
total.sim <- log(rowSums(index))
total.bts <- bts.est
total.ats <- ats.est

bts.errors <- rbind(data.frame(year=years, gear='bts', model='combined', relerror=(bts.est2-bts.sim)/bts.sim),
      data.frame(year=years, gear='bts', model='independent', relerror=(bts.est-bts.sim)/bts.sim))
ats.errors <- rbind(data.frame(year=years, gear='ats', model='combined', relerror=(ats.est2-ats.sim)/ats.sim),
      data.frame(year=years, gear='ats', model='independent', relerror=(ats.est-ats.sim)/ats.sim))

total.errors <-
  rbind(data.frame(year=years, gear='total', model='ats', relerror=(total.ats-total.sim)/total.sim),
        data.frame(year=years, gear='total', model='bts', relerror=(total.bts-total.sim)/total.sim))

all.errors <- rbind(bts.errors, ats.errors)
ggplot(all.errors, aes(year, relerror)) + geom_line()+ geom_point() +
  facet_grid(gear~model) + geom_abline(intercept=0, slope=0)

ggplot(total.errors, aes(year, relerror)) + geom_line()+ geom_point() +
  facet_grid(gear~model) + geom_abline(intercept=0, slope=0)


