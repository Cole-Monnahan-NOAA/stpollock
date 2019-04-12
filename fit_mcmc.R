## File to run the fits to the real data
chains <- 6
options(mc.cores = chains)
source('startup.R')
model <- 'combined'


## Setup a sequence of models with increasing complexity to explore when it
## starts to break

## Start with simplest model possible and do w/ and w/o the combined
## likelihood
for(combinedoff in c(TRUE,FALSE)){
control <- list(seed=121, beta2temporal=FALSE, filterdata=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                beta1temporal=FALSE, fixlambda=12, combinedoff=combinedoff)
savedir <- paste0(getwd(), '/mcmc_NS_combined',combinedoff, '_nolambda')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmc_fit.RDS'))
ind <- calculate.index.mcmc(Obj, fit)
plot.index.mcmc(ind, savedir)
plot.slow.mcmc(fit, savedir)
}

## Now add in the two catchabilities
## Start with simplest model possible and do w/ and w/o the combined
## likelihood
for(combinedoff in c(TRUE,FALSE)){
control <- list(seed=121, beta2temporal=FALSE, filterdata=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                beta1temporal=FALSE, fixlambda=0, combinedoff=combinedoff)
savedir <- paste0(getwd(), '/mcmc_NS_combined',combinedoff, '_yeslambda')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmc_fit.RDS'))
ind <- calculate.index.mcmc(Obj, fit)
plot.index.mcmc(ind, savedir)
plot.slow.mcmc(fit, savedir)
}



control <- list(seed=121, beta2temporal=FALSE, filterdata=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                beta1temporal=FALSE, fixlambda=12, combinedoff=FALSE)
savedir <- paste0(getwd(), '/mcmc_combined_nolambda_nospace')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
ind <- calculate.index.mcmc(Obj, fit)
plot.index.mcmc(ind, savedir)
plot.slow.mcmc(fit, savedir)


## This is a simplified version where the second predictor is constant in
## time but varies spatially
space <- 'NS'
control <- list(seed=121, beta2temporal=FALSE,
                beta1temporal=FALSE, temporal=2)
savedir <- paste0(getwd(), '/fit_reduced_', model, "_", space, '_', n_x)
source("prepare_inputs.R")
## Run a single iteration to optimize random effects
Obj$fn(Obj$par)
init.fn <- function() Obj$env$last.par
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=1000, open_progress=FALSE,
               init=init.fn,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir, '/fit.RDS'))
launch_shinystan(fit)

mon <- monitor(fit, print=FALSE)
pars <- names(sort(log10(mon[,'n_eff']))[1:10])
pairs(fit, pars=pars)
