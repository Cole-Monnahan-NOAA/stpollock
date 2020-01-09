

### Test sensitivity to the assumed 'inflated' zeroes for the AT
### data sets corresponding to inshore missing spatial coverage.

chains <- 5
options(mc.cores = chains)
dir.create('sensitivities/inflated0', FALSE)
td <- 13
ad <- .9
iter <- 700
warmup <- 300

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
  ##settings$FieldConfig[1:2,1:2] <- 0 # turn off spatiotemporal?
  savedir <- paste0(getwd(), '/sensitivities/inflated0/atfit_',
                    n_x,"_", zeroes.case, "/")
  fit <- fit_model(settings=settings, Lat_i=dat$Lat,
                   Lon_i=dat$Lon, t_i=dat$Year, working_dir=savedir,
                   b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2,
                   knot_method="grid", test_fit=FALSE)
  saveRDS(fit, file=paste0(savedir,'fit.RDS'))
  plot(fit, working_dir=savedir, check_residuals=FALSE)
  fits[[zeroes.case]] <- fit
}
indices <- lapply(names(fits), function(i) {
  x=fits[[i]]$parameter_estimates$SD
  ind=names(x$value) =='ln_Index_cyl'
  data.frame(case=i, est=x$value[ind], se=x$sd[ind])
}) %>% do.call(rbind,.) %>% group_by(case) %>%
  mutate(lwr=est-1.96*se, upr=est+1.96*se, year=2007:2018)
saveRDS(indices, file='results/inflated0_mle.RDS')

g <- ggplot(indices, aes(year, est, ymin=lwr, ymax=upr, fill=case, color=case)) +
  geom_ribbon(alpha=.3) + geom_line(lwd=2) + ylab("log index") + theme_bw()
ggsave('plots/sensitivity_inflated0_mle.png', width=7, height=6)
## based on these results run two MCMC sensitivities on the
## combined model (1) the basecase of a buffer and small catches
## and (2) none

#### Now run MCMC on the combined model for the two important cases
n_x <- 200
control <- list(model='combined', n_x=200,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
### Don't actually need to run this since pulled from other
### sensitivity run
## zeroes.case <- 'basecase'
## savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_',
##                   n_x,'_',zeroes.case, "/")
## source("prepare_inputs.R")
## Opt <- fit_tmb(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
##                 upper=TmbList$Upper,   savedir=savedir,
##                 newtonsteps=0, control=list(trace=1))
## results1 <- process.results(Opt, Obj, Inputs, model, space, savedir)
## plot.vastfit(results1, plotmaps=TRUE)
zeroes.case <- 'none'
control$zeroes.case <- zeroes.case
savedir <- paste0(getwd(), '/sensitivities/inflated0/senfit_',
                  n_x,'_', zeroes.case, '/')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=12512,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

x1 <- readRDS('sensitivities/kappascale/senfit_kappascale_1_combined/results.mcmc.RDS')
x2 <- readRDS('sensitivities/inflated0/senfit_200_none/results.mcmc.RDS')
out1 <- rbind(cbind(x1$index.strata, case='Inflated'),
             cbind(x2$index.strata, case='None'))
levels(out1$stratum) <- c('<0.5 m', '0.5-16 m', '>16 m')
out2 <- rbind(cbind(x1$index.gear, case='Inflated'),
             cbind(x2$index.gear, case='None'))
levels(out2$gear) <- c('Acoustic', 'Bottom Trawl', 'Total')
saveRDS(list(out1,out2), file='results/inflated0.RDS')


g1 <-  out1 %>%  ggplot(aes(year, est, fill=case, color=case,
             ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.5)+ ## geom_line(lwd=1.5)+
  facet_wrap('stratum', ncol=1, scales='free') +
  ylab('log index')+theme_bw() + theme(legend.position='none')
## ggsave('plots/sensitivity_inflated0.png', g1, width=7, height=6)
g2 <-  out2 %>%  ggplot(aes(year, est, fill=case, color=case,
             ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.5)+ ## geom_line(lwd=1.5)+
  facet_wrap('gear', ncol=1, scales='free') +
  ylab('log index')+theme_bw() + labs(y=NULL)
## ggsave('plots/sensitivity_inflated0.png', g1, width=7, height=6)

p <- plot_grid(g1, g2, rel_widths=c(1,1.3))
ggsave('plots/sensitivity_inflated0.png', p, width=7, height=6)
