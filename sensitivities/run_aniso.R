## A series of sensitivity analyses to run
chains <- 6
options(mc.cores = chains)
dir.create('sensitivities/anisofits')
td <- 16
ad <- .9
iter <- 800
warmup <- 300


### The effect of assuming anisotropic parameters. Run it fixed
### at informative values, and with assumed isotropic.
for(model in c('bts', 'ats', 'combined')){
  for(H_informative in c(TRUE, FALSE)){
    control <- list(n_x=200, H_informative=H_informative,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                model=model, make_plots=FALSE)
    if(model=='combined')
      control[c('n_eps1', 'n_eps2', 'n_omega1', 'n_omega2')] <- "IID"
    savedir <- paste0(getwd(), '/sensitivities/anisofits/senfit_aniso_informative_', H_informative,'_', model)
    source("prepare_inputs.R")
    fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                   iter=iter, open_progress=FALSE, warmup=warmup,
                   init=prior.fn, thin=1, seed=1,
                   control=list(max_treedepth=td, adapt_delta=ad))
    saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
    plot.mcmc(Obj, savedir, fit)
    res <- readRDS(paste0(savedir, '/results.mcmc.RDS'))
    res$aniso.informative <- H_informative
    res$model <- model
    res$R1_in  <- res$R2_in  <- NULL
    saveRDS(res, file=file.path(savedir, 'res.RDS'))
  }
}



results.list <- lapply(list.files('sensitivities/anisofits', full.names=TRUE),
                                 function(x) readRDS(file.path(x, 'res.RDS')))

## The indices for the independent models
out1 <- do.call(rbind, lapply(results.list, function(x)
  data.frame(model=x$model, aniso.informative=x$aniso.informative, x$index.gear))) %>%
  mutate(aniso.informative=factor(aniso.informative))
g1 <- out1 %>% filter(model !='combined') %>%
  ggplot(aes(year, est, fill=aniso.informative, color=aniso.informative, group=aniso.informative, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.5) + #geom_line(lwd=1.5)+
  facet_wrap('model', ncol=1) + ylab('log index')+ theme_bw()
ggsave('plots/sensitivity_aniso.informative_independent.png', g1, width=7, height=6)


## Look at strata in the combined model
out2 <- do.call(rbind, lapply(results.list, function(x)
  data.frame(model=x$model, aniso.informative=x$aniso.informative, x$index.strata))) %>%
  mutate(aniso.informative=factor(aniso.informative))
g2 <- out2 %>% filter(model =='combined') %>% ggplot(aes(year, est, fill=aniso.informative,
             color=aniso.informative, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.5)+ ## geom_line(lwd=1.5)+
  facet_wrap('stratum', ncol=1, scales='free') + ylab('log index')+theme_bw()
ggsave('plots/sensitivity_aniso_combined.png', g2, width=7, height=6)

saveRDS(list(out1, out2), file='results/aniso.RDS')


## library(cowplot)
## g <- plot_grid(g1,g2, nrow=2)
## save_plot('plots/sensitivity_aniso.png', g, base_width=6)
