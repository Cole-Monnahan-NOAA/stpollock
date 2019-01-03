## File to run the fits to the real data
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)


n_x <- 1000 # number of knots
for(m in 1:3){
for(s in 1:2){
model <- c('ats', 'bts', 'combined')[m]
space <- c('NS', 'S', 'ST')[s]
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
}
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
