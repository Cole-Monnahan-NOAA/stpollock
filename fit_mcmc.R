## File to run the fits to the real data
chains <- 6
options(mc.cores = chains)
source('startup.R')
td <- 15
ad <- .8
iter <- 800
warmup <- 400


## Base case for paper: combined
control <- list(model='combined', n_x=200,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=TRUE)
savedir <- paste0(getwd(), '/mcmcfit_200_combined')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
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
               init=prior.fn, seed=1245123,
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
               init=prior.fn, seed=653,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)





