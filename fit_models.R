## File to run the fits to the real data
rm(list=ls())
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)

## Test combined spatial model
n_x <- 150
model <- 'combined'; space <- 'S'
savedir <- paste0(getwd(), '/fit_', model, "_", space,  "_", n_x)
source("prepare_inputs.R")
Obj$fn()
## [1] 64255.1
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=T, loopnum=3,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)

report <- Obj$report()
report$beta1_tc
report$beta2_tc
test <- ddply(Data_Geostat, .(Gear, Year), summarize, pct0=mean(Catch_KG==0))
ggplot(test, aes(Year, pct0, group=Gear)) + geom_line()
(t(results$ParHatList$beta2_ft))



## Fit all versions of model
n_x <- 200 # number of knots
for(m in 3){
for(s in 1:3){
model <- c('ats', 'bts', 'combined')[m]
space <- c('NS', 'S', 'ST')[s]
savedir <- paste0(getwd(), '/fit_', model, "_", space, '_', n_x)
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
}
}

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
n_x <- 200
model <- 'combined'
space <- 'S'
source("prepare_inputs.R")
mcmc <- tmbstan(obj=Obj, iter=1000, chains=7, lower=TmbList$Lower,
                init='last.par.best',
                upper=TmbList$Upper, control=list(max_treedepth=12),
                open_progress=FALSE)
saveRDS(mcmc, file='mcmc.RDS')
launch_shinystan(mcmc)
