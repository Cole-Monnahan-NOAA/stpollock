### Test sensitivity to the assumed 'inflated' zeroes for the AT
### data sets corresponding to inshore missing spatial coverage.

chains <- 6
options(mc.cores = chains)
dir.create('sensitivities/inflated0', FALSE)
td <- 13
ad <- .8
iter <- 800
warmup <- 300
n_x <- 40

### These zeroes are tacked on in the load_data.R script from
### file 'ats.zeroes.RDS' so if I switch out that file and fit
### the model I can run different versions and compare them

## Same controls for all test caes
control <- list(model='combined', n_x=n_x,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=TRUE)


## base case
control$zeroes.case <- 'basecase'
savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_', n_x,'_basecase')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## sensitivity where there are more inflated zeroes right up to
## the end of the AT transects
control$zeroes.case <- 'sensitivity'
savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_', n_x,'_sensitivity')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


