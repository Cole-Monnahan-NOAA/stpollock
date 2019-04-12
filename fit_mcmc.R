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
ind <- calculate.index.mcmc(Obj, fit)
launch_shinystan(fit)
pairs(fit)


## This is a simplified version where the second predictor is constant in
## time but varies spatially
space <- 'NS'
control <- list(seed=121, beta2temporal=FALSE,
                beta1temporal=FALSE, temporal=2)
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
