
### This file runs MLE versions of the model (combined doesn't work)
source("startup.R")

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
                kappaoff=12, temporal=2, fixlambda=2, aniso=TRUE,
                make_plots=TRUE)
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



