## File to run the fits to the real data
chains <- 6
options(mc.cores = chains)
source('startup.R')
td <- 15
ad <- .9
iter <- 800
warmup <- 200


## Base case for paper: combined w/ constant catchability in P1
control <- list(model='combined', n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=TRUE)
savedir <- paste0(getwd(), '/mcmcfit_100_combined')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Base case for paper: ATS
control <- list(model='ats', n_x=200,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                make_plots=TRUE)
savedir <- paste0(getwd(), '/mcmcfit_200_ats')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Base case for paper: BTS
control <- list(model='bts', n_x=200,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                make_plots=TRUE)
savedir <- paste0(getwd(), '/mcmcfit_200_bts')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)



### lambda tests
## only lambda2
control <- list(model='combined', n_x=100, fixlambda=1,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmcfit_100_combined_lambda2')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)
## only lambda1
control <- list(model='combined', n_x=100, fixlambda=2,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmcfit_100_combined_lambda1')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


