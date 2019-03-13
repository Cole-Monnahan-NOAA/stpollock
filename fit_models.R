 ## File to run the fits to the real data
rm(list=ls())
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)
indices.to.correct <- c('ColeIndex_cy', 'ln_ColeIndex_cy', 'Index_cyl', 'ln_Index_cyl')


## Fit all versions of model
n_x <- 200 # number of knots
for(m in 3){
for(s in 2){
model <- c('ats', 'bts', 'combined')[m]
space <- c('NS', 'S', 'ST')[s]
savedir <- paste0(getwd(), '/fit_', model, "_", space, '_', n_x)
options(warn=2)
options(warn=1)
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=0, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
}
}


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
