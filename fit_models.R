
source("startup.R")

## Use the posterior medians for initial values of fixed effects
fit <- readRDS("mcmc_basecase_100/mcmcfit.RDS")
pars.fixed <- apply(as.data.frame(fit)[,-Obj$env$random], 2, median)
pars.all <- apply(as.data.frame(fit), 2, median)

## Base case for paper: combined model
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='combined',
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_100')
source("prepare_inputs.R")
Obj$env$last.par <- pars.all[-length(pars.all)]
Obj$par <- pars.fixed[-length(pars.fixed)]
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Repeat with just the BTS
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='bts',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_100_bts')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Repeat with just the ATS
control$model <- 'ats'
savedir <- paste0(getwd(), '/fit_basecase_100_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Base case for paper w/ finescale on
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=TRUE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_finscale')
source("prepare_inputs.R")
Obj$par <- pars.fixed[-length(pars.fixed)]
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)


### Fit some resolution tests on the two independent model since they're
### stable
library(snowfall)
sfInit(parallel=TRUE, cpus=12)
inputs <- expand.grid(n_x=2^(4:10), fs=c(FALSE, TRUE), model=c('ats', 'bts'))
sfExport('inputs')
results.list <- sfLapply(1:nrow(inputs), function(ii){
  fs <<- inputs$fs[ii]
  n_x <- inputs$n_x[ii]
  model <- inputs$model[ii]
  source("startup.R")
  control <<- list(seed=121, beta2temporal=TRUE,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE,
                kappaoff=0, temporal=2, fixlambda=12, make_plots=FALSE)
  savedir <<- paste0(getwd(), '/resolution_tests/test_', n_x, '_', model)
  if(fs) savedir <<- paste0(savedir, '_finescale')
  control$n_x <<- n_x; control$finescale <<- fs
  control$n_eps1 <<- control$n_eps2 <<- 0
  control$n_omega1 <<- control$n_omega2 <<- 0
  source("prepare_inputs.R")
  Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                  upper=TmbList$Upper,   savedir=savedir,
                  newtonsteps=0, control=list(trace=10))
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  results$n_x <- n_x; results$finescale <- fs
  return(results)
})

