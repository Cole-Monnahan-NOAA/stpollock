
### Sensitivities to the structure of catchability coefficients lambda. The
### base case is a single value for lambda1. Sensitivities are (1)
### time-varying lambda (with tight prior) and (2) both lambda1 and lambda2
### constant

## Base case for paper: combined w/ time-varying catchability in P1.
control <- list(model='combined', n_x=100, fixlambda=-1,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=F)
savedir <- paste0(getwd(), '/mcmcfit_100_combined_tvlambda1')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Base case for paper: combined w/ kappa1 and kappa2
control <- list(model='combined', n_x=100, fixlambda=0,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=F)
savedir <- paste0(getwd(), '/mcmcfit_100_combined_lambda12')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)
