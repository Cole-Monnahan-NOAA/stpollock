
### Test sensitivity to the assumed 'inflated' zeroes for the AT
### data sets corresponding to inshore missing spatial coverage.

chains <- 6
options(mc.cores = chains)
dir.create('sensitivities/inflated0', FALSE)
td <- 12
ad <- .8
iter <- 800
warmup <- 200
n_x <- 400

### These zeroes are tacked on in the load_data.R script from
### file 'ats.zeroes.RDS' so if I switch out that file and fit
### the model I can run different versions and compare them


### First run some MLE fits to just the AT data and are faster
### and a little easier to compare.
n_x <- 400
fits <- list()
replicateyears <- filteryears <- FALSE
efh <- 16
for(zeroes.case in c('basecase', 'sensitivity1', 'sensitivity2', 'none')){
  source('data/load_data.R')
  dat <- DF2
  dat$Gear <- factor('AT')
  dat$Catch_KG <- DF2$Catch_KG+DF3$Catch_KG
  dat$AreaSwept_km2 <- 1
  settings <- make_settings(n_x=n_x, Region="eastern_bering_sea",
                            purpose="index", bias.correct=FALSE,
                            fine_scale=TRUE)
  settings$RhoConfig[1:4] <- c(2,2,2,2)
  ## settings$FieldConfig[1:2,1:2] <- 0 # turn off spatiotemporal?
  savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_',
                    n_x,"_", zeroes.case, "/")
  fit <- fit_model(settings=settings, Lat_i=dat$Lat,
                   Lon_i=dat$Lon, t_i=dat$Year, working_dir=savedir,
                   b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2,
                   knot_method="grid", test_fit=FALSE)
  plot(fit, working_dir=savedir, check_residuals=FALSE)
  fits[[zeroes.case]] <- fit
}

indices <- lapply(names(fits), function(i) {
  x=fits[[i]]$parameter_estimates$SD
  ind=names(x$value) =='ln_Index_cyl'
  data.frame(case=i, est=x$value[ind], se=x$sd[ind])
}) %>% do.call(rbind,.) %>% group_by(case) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se, year=2007:2018)
g <- ggplot(indices, aes(year, est, ymin=lwr, ymax=upr, fill=case, color=case)) +
  geom_ribbon(alpha=.3) + geom_line(lwd=2) + ylab("log index")
ggsave('plots/sensitivity_inflated0.png', width=7, height=5)

## Same controls for all test caes
 control <- list(model='combined', n_x=n_x,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                ## n_eps1=0, n_eps2=0, n_omega2="IID", n_omega1="IID",
                make_plots=TRUE)
savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_', n_x,'_basecase')
source("prepare_inputs.R")
Opt <- fit_tmb(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results1 <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results1, plotmaps=TRUE)

control$zeroes.case <- 'sensitivity1'
savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_', n_x,'_sensitivity1')
source("prepare_inputs.R")
Opt <- fit_tmb(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results2 <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results2, plotmaps=TRUE)

control$zeroes.case <- 'sensitivity2'
savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_', n_x,'_sensitivity2')
source("prepare_inputs.R")
Opt <- fit_tmb(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=F,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results3 <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results3, plotmaps=TRUE)

control$zeroes.case <- 'none'
savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_', n_x,'_none')
source("prepare_inputs.R")
Opt <- fit_tmb(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=F,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results4 <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results4, plotmaps=TRUE)



index <- rbind(cbind(version='basecase', results1$Index),
      cbind(version='low', results2$Index),
      cbind(version='high', results3$Index))
filter(index, year==2008)

ggplot(index, aes(year, est, ymin=lwr, ymax=upr, fill=version)) + geom_ribbon(alpha=.5)



fit <- fit_model(settings=settings, Lat_i=dat$lat,
                 Lon_i=dat$lon, t_i=dat$year,
                 ## observations_LL=cbind(Lon=dat$lon, Lat=dat$lat),
                 b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2)
plot(fit)



## base case
control$zeroes.case <- 'basecase'
savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_', n_x,'_basecase')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## sensitivity where there are more inflated zeroes right up to
## the end of the AT transects
control$zeroes.case <- 'sensitivity1'
savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_', n_x,'_sensitivity1')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


