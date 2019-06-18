### This files runs the simulation testing component of the analysis.

rm(list=ls())
source("startup.R")
## We condition the OM on the fitted base case model
load('fit_basecase_100/Save.RData')
load('fit_basecase_100/Record.RData')
names(Save)
names(Record)
## Rebuild the Obj with the MLEs
mle <- Save$est$est
par.names <- as.character(Save$est$par)
control <- Record$control
savedir <- paste0(getwd(), '/simulations/sim_basecase_100')
source("prepare_inputs.R")
Obj$par <- mle
Obj$fn(mle) -  Save$Opt$objective ## check MLE fit
Obj$gr(mle) ## check that the gradients are 0
pars.all <- Obj$env$last.par ## the full parameter vector

## Build a simulated OM by changing the betas
beta1_ft <- (pars.all[grep('beta1_ft', names(pars.all))])
beta1_ft <- matrix(beta1_ft, ncol=3, byrow=TRUE)
par(mfrow=c(2,2))
matplot(beta1_ft)
beta1 <- cbind(c(seq(0,-1.5, len=6), seq(-1.25, 3.5, len=6)),
               c(seq(.5,-1, len=6), seq(-1.25, 2.5, len=6)),
               c(seq(.5,-.5, len=6), seq(-.75,2, len=6)))
matplot(beta1)
pars.all[grep('beta1_ft', names(pars.all))] <- as.vector(t(beta1))
beta2_ft <- (pars.all[grep('beta2_ft', names(pars.all))])
beta2_ft <- matrix(beta2_ft, ncol=3, byrow=TRUE)
matplot(beta2_ft)
beta2 <- cbind(seq(0,1, len=12), seq(0, -.5, len=12), seq(.5, -.5, len=12))
matplot(beta2)
pars.all[grep('beta2_ft', names(pars.all))] <- as.vector(t(beta2))
## Put data into the ATS for missing years so it's sampled in the
## simulation
DF1 <- subset(Data_Geostat, Gear=='Trawl')
DF2 <- subset(Data_Geostat, Gear=='Acoustic_3-16')
DF3 <- subset(Data_Geostat, Gear=='Acoustic_16-surface')
tmp2 <- subset(Data_Geostat, Gear == 'Acoustic_3-16' & Year %in% c(2010, 2012, 2014, 2016))
tmp2$Year <- tmp2$Year+1
DF2 <- rbind(DF2, tmp2)
tmp3 <- subset(Data_Geostat, Gear == 'Acoustic_16-surface' & Year %in% c(2010, 2012, 2014, 2016))
tmp3$Year <- tmp3$Year+1
DF3 <- rbind(DF3, tmp3)
## Have to rebuild the object since the data changed
control$simdata <- TRUE
source("prepare_inputs.R")
table(Data_Geostat$Gear)
## Now when I call a new vector it'll have expected values in the missing
## years and thus be sampled from them.


### Loop through replicates of the simulation
## for(i in 1:3){
i <- 1
## Simulate the sampling process.
set.seed(i)
sim <- Obj$report(pars.all)
encounters <- rbinom(length(sim$R1_i), size=1, prob=sim$R1_i)
## Carefully index to get the right SigmaObs by gear type
sigma.obs <- sim$SigmaM[,1][as.numeric(Data_Geostat$Gear)]
logcatches <- rnorm(length(sim$R2_i), log(sim$R2_i)-sigma.obs^2/2, sigma.obs)
catches <- exp(logcatches)*encounters

## Rebuild the Obj with the new simulated data
Data_Geostat$Catch_KG <- catches
DF1 <- subset(Data_Geostat, Gear=='Trawl')
DF2 <- subset(Data_Geostat, Gear=='Acoustic_3-16')
DF3 <- subset(Data_Geostat, Gear=='Acoustic_16-surface')

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
Obj$par <- mle # cheat by starting at MLE
tryCatch(Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=getsd,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1)), error='XX')
if(is.list(Opt)){
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  plot.vastfit(results, plotmaps=TRUE)
}

## Repeat with just the BTS
control$model <- 'bts'; control$make_plots <- FALSE
savedir <- paste0(getwd(), '/simulations/simfit_', i, '_bts')
source("prepare_inputs.R")
saveRDS(sim, paste0(savedir, '/simreport.RDS'))
tryCatch(Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=getsd,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1)), error='XX')
if(is.list(Opt)){
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  plot.vastfit(results, plotmaps=TRUE)
}

## Repeat with just the ATS
control$model <- 'ats'
savedir <- paste0(getwd(), '/simulations/simfit_', i, '_ats')
source("prepare_inputs.R")
saveRDS(sim, paste0(savedir, '/simreport.RDS'))
tryCatch(Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=getsd,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1)), error='XX')
if(is.list(Opt)){
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  plot.vastfit(results, plotmaps=TRUE)
}



## Process and compare the fits to the simulated truth
index.sim <- sim$ln_ColeIndex_cy
index <- t(sim$Index_cyl[,,1])
matplot(log(index))
fitc <- readRDS
load('simulations/simfit_combined_1/Save.RData')
indexc <- t(Save$Report$Index_cyl[,,1])
load('simulations/simfit_bts_1/Save.RData')
indexb <- Save$Report$Index_cyl[,,1]
load('simulations/simfit_ats_1/Save.RData')
indexa <- Save$Report$Index_cyl[,,1]

bts.sim <- log(rowSums(index[, 1:2]))
bts.est <- log(indexb)
bts.est2 <- log(rowSums(indexc[,1:2]))
ats.sim <- log(rowSums(index[, 2:3]))
ats.est <- log(indexa)
ats.est2 <- log(rowSums(indexc[,2:3]))

bts.errors <- rbind(data.frame(year=years, gear='bts', model='combined', relerror=(bts.est2-bts.sim)/bts.sim),
      data.frame(year=years, gear='bts', model='independent', relerror=(bts.est-bts.sim)/bts.sim))
ats.errors <- rbind(data.frame(year=years, gear='ats', model='combined', relerror=(ats.est2-ats.sim)/ats.sim),
      data.frame(year=years, gear='ats', model='independent', relerror=(ats.est-ats.sim)/ats.sim))

all.errors <- rbind(bts.errors, ats.errors)
ggplot(all.errors, aes(year, relerror)) + geom_line()+ geom_point() +
  facet_grid(gear~model) + geom_abline(intercept=0, slope=0)


