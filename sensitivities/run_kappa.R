
## A series of sensitivity analyses to run
chains <- 6
options(mc.cores = chains)
td <- 16
ad <- .9
iter <- 1000
warmup <- 300
dir.create('sensitivities/kappascale')

## I originally ran this on the independent ones but I don't
## think it's really necessary and just takes time so just run
## the combined model

### Note that the kappascale=.5 is REALLY slow because of big
### tree depths. That's why td is 16 here

### The effect of fixing logkappa. Run the models with quarter and
### quadruple the spatial range used
model <- 'combined'
for(kappascale in c(.5,1,2)){
  control <- list(n_x=200, n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                  model=model, kappascale=kappascale)
  if(model=='combined')
    control[c('n_eps1', 'n_eps2', 'n_omega1', 'n_omega2')] <- "IID"
  savedir <- paste0(getwd(), '/sensitivities/kappascale/senfit_kappascale_', kappascale,'_', model)
  source("prepare_inputs.R")
  fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                 iter=iter, open_progress=FALSE, warmup=warmup,
                 init=prior.fn, thin=1, seed=1,
                 control=list(max_treedepth=td, adapt_delta=ad))
  saveRDS(object = fit, file=paste0(savedir,'/senfit.RDS'))
  plot.mcmc(Obj, savedir, fit)
  ## Only save the pieces needed
  index.strata <- readRDS(file.path(savedir, 'results.mcmc.RDS'))$index.strata
  res <- list(index.strata=index.strata,
              kappascale=kappascale, model=model)
  saveRDS(res, file=file.path(savedir, 'res.RDS'))
}



results.list <-
  lapply(list.files('sensitivities/kappascale', full.names=TRUE,
                    pattern='kappascale'),
         function(x) readRDS(file.path(x, 'res.RDS')))
## Look at strata in the combined model
out <- do.call(rbind, lapply(results.list, function(x)
  data.frame(model=x$model, kappascale=x$kappascale, x$index.strata))) %>%
  mutate(kappascale=factor(kappascale))

levels(out$stratum) <- c('<0.5 m', '0.5-16 m', '>16 m')

g <- out %>% filter(model =='combined') %>%
  ggplot(aes(x=year, y=est, fill=kappascale,   group=kappascale, ymin=lwr, ymax=upr)) +
  geom_ribbon(alpha=.3) + geom_line(aes(color=kappascale), alpha=.6, lwd=1, show.legend=FALSE)+
  facet_wrap('stratum', ncol=1, scales='free') +
  ylab('log index')+theme_bw() + labs(fill='kappa\nmultiplier')
ggsave('plots/sensitivity_kappascale.png', g, width=7, height=6)

saveRDS(out, file='results/kappascale.RDS')



