## A series of sensitivity analyses to run
chains <- 1
options(mc.cores = chains)
setwd('..')
source('startup.R')

### The effect of fixing logkappa. Run the models with half and double the
### spatial range used (50km and 200km)
control <- list(beta2temporal=TRUE, n_x=100,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE,
                kappaoff=12, temporal=2)
for(model in c('bts', 'ats', 'combined')){
  for(kappascale in c(.5,1,2)){
    control$model <- model
    control$kappascale <- kappascale
    savedir <- paste0(getwd(), '/sensitivities/kappascale_', kappascale,'_', model)
    source("prepare_inputs.R")
    fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                   iter=800, open_progress=FALSE, warmup=200,
                   init='last.par.best', thin=1,
                   control=list(max_treedepth=12))
    saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
    plot.mcmc(Obj, savedir, fit)
  }
}


