## File to run the fits to the real data
chains <- 8
source('startup.R')

## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)

## Fit 2 versions of combined model with tmbstan
model <- c('ats', 'bts', 'combined')[3]
#space <- c('NS', 'S', 'ST')[3]
n_x <- 80
space <- 'NS'
control <- list(seed=121, temporal=2, beta2temporal=FALSE)
savedir <- paste0(getwd(), '/fit_full_', model, "_", space, '_', n_x)
source("prepare_inputs.R")
## Run a single iteration to optimize random effects
Obj$fn(Obj$par)
init.fn <- function() Obj$env$last.par
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=1000, open_progress=FALSE,
               init=init.fn,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/fit.RDS'))
launch_shinystan(fit)


## This is a simplified version where the second predictor is constant in
## time but varies spatially
space <- 'S'
control <- list(seed=121, n_eps1=0, n_eps2=0, beta2temporal=FALSE, temporal=2)
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





## Test bias adjustment for index. Run this once then again with a slightly different n_x and change
## prepare_inputs to have epsilon rho = 0. Thus there's w & w/o bias correction on a model with and without a
## smoother on ST effects
model <- 'combined'
space <- 'ST'
n_x <- 200
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
