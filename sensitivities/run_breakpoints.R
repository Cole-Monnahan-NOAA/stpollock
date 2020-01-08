

## Look at how the ATS changes with different breakpoints for the lower and
## upper bounds (EFH). We have two cases: [.5,16] (Stan's EFH
## estimate) and [.5,3] which implies no vertical schooling so
## the EFH is the physical height of the gear.

## Note that there are two input data files created specifically
## for these. It cannot be changed from 3 and 16. The control$efh
## changes which file is read in and used.

chains <- 5
options(mc.cores = chains)
td <- 13
ad <- .9
iter <- 700
warmup <- 200
dir.create('sensitivities/breakpoints', showWarnings=FALSE)


control <- list(model='combined', n_x=100,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)

## Case 1: base case from paper: take this from sensitivity tests
## control$efh <- 16
## savedir <- paste0(getwd(), '/sensitivities/breakpoints/senfit_efh16')
## source("prepare_inputs.R")
## fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
##                iter=iter, open_progress=FALSE, warmup=warmup,
##                init=prior.fn, seed=85234,
##                control=list(max_treedepth=td, adapt_delta=ad))
## saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
## plot.mcmc(Obj, savedir, fit)

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


x1 <- readRDS('sensitivities/kappascale/senfit_kappascale_1_combined/results.mcmc.RDS')
x2 <- readRDS('sensitivities/breakpoints/senfit_efh3/results.mcmc.RDS')
index.gear <- rbind(cbind(x1$index.gear, EFH=factor(16)), cbind(x2$index.gear, EFH=factor(3)))
index.strata <- rbind(cbind(x1$index.strata, EFH=factor(16)), cbind(x2$index.strata, EFH=factor(3)))
levels(index.gear$gear) <- c('Acoustic', 'Bottom Trawl', 'Total')
levels(index.strata$stratum) <- c('<0.5 m', '0.5-EFH', '>EFH')
index <- list(index.gear, index.strata)
saveRDS(index, file='results/breakpoints.RDS')

g1 <- ggplot(index.strata, aes(year, est, ymin=lwr, ymax=upr, color=EFH, fill=EFH)) +
  geom_line(lwd=2)+
  geom_ribbon(alpha=1/3) + ylab('log index') +
  facet_wrap('stratum', ncol=1) + theme_bw() +
  theme(legend.position = 'none') + ggtitle('By depth layer')
g2 <- ggplot(index.gear, aes(year, est, ymin=lwr, ymax=upr, color=EFH, fill=EFH)) +
  geom_line(lwd=2)+
  geom_ribbon(alpha=1/3) + ylab('') + ggtitle('By gear type')+
  facet_wrap('gear', ncol=1) + theme_bw()
## combine those two above into one
library(cowplot)
g <- cowplot::plot_grid(g1, g2, ncol=2, rel_widths=c(1,1.2))
ggsave('plots/sensitivity_breakpoints.png', g, width=9, height=7)

