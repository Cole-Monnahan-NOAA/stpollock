
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
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='ats',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
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





control <- list(seed=121, beta2temporal=TRUE, n_x=51,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE,
                kappaoff=12, temporal=0, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_base_ST')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=10))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Base + smoothing even the no missing years
control <- list(seed=121, beta2temporal=TRUE, n_x=51,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_smoother_ST')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=10))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Base + smoothing w/ missing years
control <- list(seed=121, beta2temporal=TRUE, n_x=51,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_smootherfull_ST')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=10))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)


## Base model + smoothing on years but  w/ missing years + catchability
## estimated + more knots but no finescale
control <- list(seed=121, beta2temporal=TRUE, n_x=200,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/fit_smootherfull_lambda1_hires')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=10))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)


## Test bias adjustment for index. Run this once then again with a slightly different n_x and change
## prepare_inputs to have epsilon rho = 0. Thus there's w & w/o bias correction on a model with and without a
## smoother on ST effects
savedir <- paste0(getwd(), '/bias_', model, "_", space, '_', n_x)
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=5,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=10))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
## Repeat with bias adjust turned on for index
savedir <- paste0(getwd(), '/bias_cor_', model, "_", space, '_', n_x)
source("prepare_inputs.R")
Obj$par <- Opt$par
indices.to.correct <- c('ColeIndex_cy', 'ln_ColeIndex_cy', 'Index_cyl', 'ln_Index_cyl')
Opt2 <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10),
                bias.correct=TRUE,
                bias.correct.control=list(vars_to_correct=indices.to.correct))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt2, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
cbind(Opt$time_for_run, Opt$time_for_sdreport )
cbind(Opt2$time_for_run, Opt2$time_for_sdreport )


## Test increasing resolution
for(n_x in 2^(11:12)){
  space <- "S"; model <- 'combined'
  savedir <- paste0(getwd(), '/knots_combined_S_',n_x)
  source("prepare_inputs.R")
  Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                  upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
  ## TMBhelper::Check_Identifiable(Obj)
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
}


### Some MCMC tests
library(tmbstan)
library(shinystan)
options(mc.cores = 7)
## n_x <- 200
## model <- 'combined'
## space <- 'S'
## source("prepare_inputs.R")
lwr <- TmbList$Lower
lwr[grep('L_', names(lwr))] <- 0
mcmc <- tmbstan(obj=Obj, iter=1000, chains=7,
                init='par', upper= TmbList$Upper,
                lower=lwr,
                 control=list(max_treedepth=8, adapt_delta=.9),
                open_progress=FALSE)
saveRDS(mcmc, file='mcmc.RDS')
launch_shinystan(mcmc)
