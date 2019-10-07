
## File to run the sensitivities for spatial configurations
chains <- 6
options(mc.cores = chains)
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
               init=prior.fn, seed=62,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Spatial only
control <- list(model='combined', n_x=100,
                n_eps1=0, n_eps2=0, n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/spatialconfig/senfit_S')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=5185312,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Spatiotemporal
control <- list(model='combined', n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/spatialconfig/senfit_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=87153,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)



results.list <- lapply(list.files('sensitivities/spatialconfig', full.names=TRUE),
                                 function(x) readRDS(file.path(x, 'results.mcmc.RDS')))
## The indices for the independent models
out <- do.call(rbind, lapply(results.list, function(x)
  cbind(x$index.strata, spatialconfig=x$index.gear$space[1])))  %>%
  mutate(Configuration=factor(spatialconfig, levels=c("NS", "S",
                                 "ST"), labels=c("No Space",
                                 "Spatial", 'Spatiotemporal')))
saveRDS(out, file='results/spatialconfig.RDS')
g1 <- ggplot(out, aes(year, est, fill=Configuration, color=Configuration, group=Configuration, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=1/3) +# geom_line(lwd=1.5)+
  facet_wrap('stratum', ncol=1, scales='free_y') + ylab('log index')+ theme_bw()
ggsave('plots/sensitivity_spatialconfig.png', g1, width=7, height=6)





