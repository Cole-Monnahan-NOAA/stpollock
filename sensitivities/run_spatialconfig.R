
## File to run the sensitivities for spatial configurations
chains <- 6
options(mc.cores = chains)
setwd('..')
source('startup.R')
dir.create('sensitivities/spatialconfig')
td <- 15
ad <- .9
iter <- 800
warmup <- 400

## Combined with no spatial
control <- list(model='combined', n_x=100,
                n_eps1=0, n_eps2=0, n_omega2=0, n_omega1=0,
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/spatialconfig/senfit_NS')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

control <- list(model='combined', n_x=100,
                n_eps1=0, n_eps2=0, n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/spatialconfig/senfit_S')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

control <- list(model='combined', n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/spatialconfig/senfit_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)



results.list <- lapply(list.files('sensitivities/spatialconfig', full.names=TRUE),
                                 function(x) readRDS(file.path(x, 'results.mcmc.RDS')))

## The indices for the independent models
out <- do.call(rbind, lapply(results.list, function(x)
  cbind(x$index.strata, spatialconfig=x$index.gear$space[1])))  %>%
#  data.frame(model=x$model, spatialconfig=x$spatialconfig, x$index.gear))) %>%
  mutate(spatialconfig=factor(spatialconfig))
g1 <- ggplot(out, aes(year, est, fill=spatialconfig, color=spatialconfig, group=spatialconfig, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.5) + geom_line(lwd=1.5)+
  facet_wrap('stratum', ncol=1, scales='free_y') + ylab('log index')+ theme_bw()
ggsave('sensitivities/sensitivity_spatialconfig.png', g1, width=7, height=6)





