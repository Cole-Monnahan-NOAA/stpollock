
## Look at how the ATS changes with different breakpoints for the lower and
## upper bounds (EFH). We have two cases: [.5,16] (Stan's EFH
## estimate) and [.5,3] which implies no vertical schooling so
## the EFH is the physical height of the gear.

## Note that there are two input data files created specifically
## for these. It cannot be changed from 3 and 16. The control$efh
## changes which file is read in and used.

chains <- 6
options(mc.cores = chains)
td <- 12
ad <- .8
iter <- 800
warmup <- 200
dir.create('sensitivities/breakpoints', showWarnings=FALSE)


control <- list(model='combined', n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
## The flag "fixlambda" controls the catchability
## configuration. Kind of a bad way to set it up, see
## prepare.inputs for how it actually works.

## Case 1: base case from paper
control$efh <- 16
savedir <- paste0(getwd(), '/sensitivities/breakpoints/senfit_efh16')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=85234,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Case 2: sensitivity of 3m
control$efh <- 3
savedir <- paste0(getwd(), '/sensitivities/breakpoints/senfit_efh3')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, seed=85234,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


x1 <- readRDS('sensitivities/breakpoints/senfit_efh16/results.mcmc.RDS')$index.gear
x2 <- readRDS('sensitivities/breakpoints/senfit_efh3/results.mcmc.RDS')$index.gear
index.gear <- rbind(cbind(x1, EFH=factor(16)), cbind(x2, EFH=factor(3)))
x1 <- readRDS('sensitivities/breakpoints/senfit_efh16/results.mcmc.RDS')$index.strata
x2 <- readRDS('sensitivities/breakpoints/senfit_efh3/results.mcmc.RDS')$index.strata
index.strata <- rbind(cbind(x1, EFH=factor(16)), cbind(x2, EFH=factor(3)))
index <- list(index.gear, index.strata)

saveRDS(index, file='results/breakpoints.RDS')

g1 <- ggplot(index.strata, aes(year, est, ymin=lwr, ymax=upr, color=EFH, fill=EFH)) +
  geom_line(lwd=2)+
  geom_ribbon(alpha=1/3) + ylab('log index') +
  facet_wrap('stratum', ncol=1) + theme_bw()

g2 <- ggplot(index.gear, aes(year, est, ymin=lwr, ymax=upr, color=EFH, fill=EFH)) +
  geom_line(lwd=2)+
  geom_ribbon(alpha=1/3) + ylab('log index') +
  facet_wrap('gear', ncol=1) + theme_bw()

## combine those two above into one
library(cowplot)
g <- cowplot::plot_grid(g1, g2, ncol=2)
ggsave('plots/sensitivities_breakpoints.png', g, width=9, height=7)

