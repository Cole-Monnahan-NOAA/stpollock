## File to run the fits to the real data
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)


n_x <- 300 # number of knots
for(m in 1:3){
for(s in 2){
model <- c('ats', 'bts', 'combined')[m]
space <- c('NS', 'S', 'ST')[s]
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, model, space, savedir)
plot.vastfit(results)
}
}


### Some MCMC tests
library(tmbstan)
options(mc.cores = 4)
n_x <- 400
model <- 'combined'
space <- 'S'
source("prepare_inputs.R")
mcmc <- tmbstan(obj=Obj, iter=1000, chains=4, lower=TmbList$Lower, upper=TmbList$Upper)
