## A series of sensitivity analyses to run
chains <- 6
options(mc.cores = chains)
td <- 15
ad <- .9
iter <- 800
warmup <- 400
dir.create('sensitivities/kappascalefits')

## I originally ran this on the independent ones but I don't
## think it's really necessary and just takes time so just run
## the combined model


### The effect of fixing logkappa. Run the models with half and
### double the spatial range used
for(model in c('bts', 'ats', 'combined')[3]){
  for(kappascale in c(.25,1,4)){
    control <- list(n_x=100, n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                    model=model, kappascale=kappascale)
    if(model=='combined')
      control[c('n_eps1', 'n_eps2', 'n_omega1', 'n_omega2')] <- "IID"
    savedir <- paste0(getwd(), '/sensitivities/kappascalefits/mcmcfit_kappascale_', kappascale,'_', model)
    source("prepare_inputs.R")
    fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                   iter=iter, open_progress=FALSE, warmup=warmup,
                   init=prior.fn, thin=1, seed=1,
                   control=list(max_treedepth=td, adapt_delta=ad))
    saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
    res <- get.results.mcmc(Obj, fit)
    res$kappascale <- kappascale
    res$logkappainput1 <- logkappainput1
    res$logkappainput2 <- logkappainput2
    res$model <- model
    ## These are big and not needed
    res$R1_in  <- res$R2_in  <- NULL
    saveRDS(res, file=file.path(savedir, 'res.RDS'))
  }
}

results.list <- lapply(list.files('sensitivities/kappascalefits', full.names=TRUE,
                                  pattern='kappascale'), function(x)
                       readRDS(file.path(x, 'res.RDS')))

out1 <- NULL
## ## The indices for the independent models
## out1 <- do.call(rbind, lapply(results.list, function(x)
##   data.frame(model=x$model, kappascale=x$kappascale, x$index.gear))) %>%
##   mutate(kappascale=factor(kappascale))
## g1 <- out1 %>% filter(model !='combined') %>%
##   ggplot(aes(year, est, fill=kappascale, color=kappascale, group=kappascale, ymin=lwr, ymax=upr)) +
##   geom_ribbon(alpha=.3) + geom_line(lwd=1.5)+
##   facet_wrap('model', ncol=1, scales='free') + ylab('log index')+ theme_bw()
## ggsave('plots/sensitivity_kappascale_independent.png', g1, width=7, height=6)


## Look at strata in the combined model
out2 <- do.call(rbind, lapply(results.list, function(x)
  data.frame(model=x$model, kappascale=x$kappascale, x$index.strata))) %>%
  mutate(kappascale=factor(kappascale))
g2 <- out2 %>% filter(model =='combined') %>%
  ggplot(aes(year, est, fill=kappascale, color=kappascale, group=kappascale, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.3) + geom_line(lwd=1.5)+
  facet_wrap('stratum', ncol=1, scales='free') +
  ylab('log index')+theme_bw()
ggsave('plots/sensitivity_kappascale_combined.png', g2, width=7, height=6)

saveRDS(list(out1,out2), file='results/kappa.RDS')

## library(cowplot)
## g <- plot_grid(g1,g2, nrow=2)
## save_plot('plots/sensitivity_kappascale.png', g, base_width=6)


