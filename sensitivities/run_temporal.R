## Explore the difference between AR1 and random walk for the
## combined model. Note that I have informative priors built
## into the model for the rho parameter that turn on when the par
## is active only

td <- 15
ad <- .8
chains <- 6
options(mc.cores = chains)
iter <- 800
warmup <- 300
dir.create('sensitivities/temporal', showWarnings=FALSE)

## Setup with AR1 or RW
for(tmp in c('AR1', 'RW')){
  control <- list(model='combined', n_x=100,
                  temporal=ifelse(tmp=='AR1', 4,2),
                  ## n_eps1=0, n_eps2=0, n_omega2=0, n_omega1=0,
                  n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                  make_plots=FALSE)
  savedir <- paste0(getwd(), '/sensitivities/temporal/senfit_', tmp)
  source("prepare_inputs.R")
  fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                 iter=iter, open_progress=FALSE, warmup=warmup,
                 init=prior.fn, thin=1, seed=1,
                 control=list(max_treedepth=td, adapt_delta=ad))
  saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
  plot.mcmc(Obj, savedir, fit)
}


ar1 <- readRDS('sensitivities/temporal/senfit_AR1/results.mcmc.RDS')
rw <- readRDS('sensitivities/temporal/senfit_RW/results.mcmc.RDS')
res <- rbind(cbind(temporal='AR1', ar1$index.strata),
             cbind(temporal='RW', rw$index.strata))
## Index comparison
g <- ggplot(res, aes(year, est, ymin=lwr, ymax=upr, fill=temporal)) +
  geom_ribbon(alpha=.5) +
  facet_wrap('stratum', ncol=1, scales='free') +
  ylab('log index')+theme_bw()
ggsave('plots/sensitivity_temporal.png', g, width=7, height=6)

## plot marginal posteriors of rhos vs the prior
rhos <- readRDS('sensitivities/temporal/senfit_AR1/mcmcfit.RDS') %>%
  as.data.frame() %>% select(contains('rho')) %>%
  pivot_longer(cols=1:4) %>%
  separate(name, into=c('par', 'lp', 'layer'),  sep='_') %>%
  mutate(layer=factor(layer, levels=c('f[1]', 'f[2]', 'f[3]'),
                      labels=c('<0.5m', '0.5-16m', '>16m')),
         Rho=case_when(lp=='rho1'~'rho n',
                       lp=='rho2'~'rho w'))
## Prior is rho~dbeta(2,2) so recreate that to add
xseq<- seq(-.999,.999, len=1000)
prior <- expand.grid(par=unique(rhos$par), lp=unique(rhos$lp),
                     layer=unique(rhos$layer),
                     value=xseq) %>%
  ## Divide by two b/c of Jacobian
  mutate(density=dbeta((1+value)/2, 2,2)/2)
g <- ggplot(rhos, aes(value, after_stat(density), fill=Rho)) +
  geom_histogram(alpha=.5, , position='identity') +
  geom_line(data=prior, aes(x=value, y=density, group=factor(lp), fill=NULL)) +
  facet_grid(par~.) + theme_bw() + xlim(-1,1) + labs(x='Rho')
ggsave('plots/sensitivity_temporal_rho_posterior.png', g, width=7, height=6)

