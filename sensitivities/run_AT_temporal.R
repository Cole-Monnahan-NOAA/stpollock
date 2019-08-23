## Explore the temporal cnofiguration for the AT data set

setwd('..')
source("startup.R")

## Repeat with just the ATS
control <- list(model='ats', n_x=100,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1, make_plots=TRUE)
savedir <- paste0(getwd(), '/senfit_100_ats_ar1')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

control <- list(model='ats', n_x=100,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1, make_plots=FALSE)
control$temporal <- 2 ## random walk
savedir <- paste0(getwd(), '/senfit_100_ats_rw')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Repeat with just the ATS
control <- list(model='ats', n_x=100, replicateyears=TRUE,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1, make_plots=TRUE)
savedir <- paste0(getwd(), '/senfit_ats_ar1_replicateyears')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

control <- list(model='ats', n_x=100, replicateyears=TRUE,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1, make_plots=FALSE)
control$temporal <- 2 ## random walk
savedir <- paste0(getwd(), '/senfit_ats_rw_replicateyears')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
