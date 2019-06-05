## File to run the fits to the real data
chains <- 10
options(mc.cores = chains)
source('startup.R')
model <- 'combined'


## Base case for paper
control <- list(seed=121, beta2temporal=TRUE, n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/mcmc_basecase_100')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=800, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## Base model: IID ST w/ kappas off and filtered years but no temporal
## smoothing on anything.
control <- list(seed=121, beta2temporal=TRUE, n_x=150,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE,
                kappaoff=12, temporal=0, fixlambda=12, make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmc_base_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Base model + smoothing on years but no missing years
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE,
                kappaoff=12, temporal=2, fixlambda=12, make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmc_smoother_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best', thin=1,
               control=list(max_treedepth=10))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Base model + smoothing on years but  w/ missing years
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE,
                kappaoff=12, temporal=2, fixlambda=12, make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmc_smootherfull_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## Base model + smoothing on years but  w/ missing years + catchability estimated
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmc_smootherfull_lambda1_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE, warmup=200,
               init='last.par.best', thin=1,
               control=list(max_treedepth=8))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## Base model + smoothing on years but  w/ missing years + catchability
## estimated + finescale
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE, finescale=TRUE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/mcmc_smootherfull_lambda1_finscale_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best', thin=1,
               control=list(max_treedepth=12))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## Temp code to plot indices and availability by the three scenarios above
x1 <- readRDS('mcmc_kappa1est_fixlambda1_ST/index.mcmc.RDS')
x2 <- readRDS('mcmc_kappaoff_tvlambda_ST/index.mcmc.RDS')
x3 <- readRDS('mcmc_kappa1est_tvlambda_ST/index.mcmc.RDS')
index.gear <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='By Gear', x$index.gear, scenario=x$scenario)))
index.strata <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='By Strata', x$index.strata, scenario=x$scenario)))
avail <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='Availability', x$availability, scenario=x$scenario)))
library(cowplot)
g1 <- ggplot(index.gear, aes(year, y=est, color=gear, group=gear, fill=gear)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('log abundance')
g2 <- ggplot(index.strata, aes(year, y=est, color=stratum, group=stratum, fill=stratum)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('log abundance')
g3 <- ggplot(avail, aes(year, y=est, color=gear, group=gear, fill=gear)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('Availability to gear')
g <- plot_grid(g1, g2,g3, ncol = 1)
save_plot("plots/scenario_grid.png", g, base_width=9, base_height=6,
          base_aspect_ratio = 1.3)

