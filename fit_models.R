## File to run the fits to the real data
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)


n_x <- 500 # number of knots
model <- c('ats', 'bts', 'combined')[3]
space <- c('NS', 'S', 'ST')[3]
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, model, space, savedir)
plot.vastfit(results)
