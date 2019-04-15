## File to run the fits to the real data
chains <- 8
options(mc.cores = chains)
source('startup.R')
model <- 'combined'


## Start with simplest model possible and do w/ and w/o the combined
## likelihood and different combinations of estimating catchability
for(combinedoff in c(TRUE,FALSE)){
for(fixlambda in c(0,1,2,12)){
control <- list(seed=121, beta2temporal=FALSE, filterdata=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                beta1temporal=FALSE, fixlambda=fixlambda,
                combinedoff=combinedoff )
savedir <- paste0(getwd(), '/mcmc_NS_combined',combinedoff, '_lambda',fixlambda)
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE, thin=2,
               init='last.par.best',
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmc_fit.RDS'))
plot.mcmc(Obj, savedir, fit)
}
}


## This is our simplest base case model with the subsetted data (no
## temporal aspect
control <- list(seed=121, beta2temporal=FALSE, n_x=50,
                beta1temporal=FALSE, temporal=2, n_eps2=0)
savedir <- paste0(getwd(), '/mcmc_combined_ST')
source("prepare_inputs.R")
## Run a single iteration to optimize random effects
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=14))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)
ind <- calculate.index.mcmc(Obj, fit)
plot.index.mcmc(ind, savedir)
plot.slow.mcmc(fit, savedir)
