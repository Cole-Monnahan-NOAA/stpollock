### This files runs the simulation testing component of the analysis.

rm(list=ls())
source("startup.R")
nsim <- 100
n_x <- 300 # number of knots for OM and EM (they match)
clean.dir <- function(savedir){
  ## Unlink and delete TMB objects to save space and prevent errors about
  ## too many dlls loaded.
  dyn.unload(paste0(savedir,'/VAST_v8_0_0.dll'))
  x <- paste0(savedir,'/', c('VAST_v8_0_0.dll', 'VAST_v8_0_0.o', 'Record.RData'))
  trash <- file.remove(x)
}

### Step 1: Build a base OM model from which to simulate
### data. This model has no ST effects and independent spatial
### effects among strata. We use the anisotropy from pcod in the
### OM and estimate these parameters in the EM. We also use the
### kappas from pcod and estimate these.
control <- list(beta2temporal=TRUE, n_x=n_x, model='combined',
                ##n_eps1=0, n_eps2=0, n_omega2=0, n_omega1=0,
                n_eps1=0, n_eps2=0, n_omega2="IID", n_omega1='IID',
                beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                kappaoff=0, temporal=0, fixlambda=2,
                simdata=FALSE, simulation=TRUE,
                aniso=TRUE)
savedir <- paste0(getwd(), '/simulations/OM')
source("prepare_inputs.R")
Obj.OM <- Obj
par.names <- names(Obj.OM$env$last.par)
par.truth0 <- Obj.OM$env$last.par*NA ## make sure everything is changed
dat0 <- Data_Geostat ## save a copy of original data since gets overwritten
nyrs <- length(unique(Data_Geostat$Year))
beta1.trend <- cbind(c(seq(2,-1, len=nyrs-3), -1,0,1),
               c(seq(-1, .5, len=nyrs-3), .5,.5,.5),
               seq(-1,1, len=nyrs))
beta2.trend <- cbind(seq(4,7, len=nyrs),
               seq(5, 6, len=nyrs),
               seq(4, 6, len=nyrs))
beta1.flat <- cbind(seq(1,1, len=nyrs),
               seq(-1, -1, len=nyrs),
               seq(-2,-2, len=nyrs))
beta2.flat <- cbind(seq(5,5, len=nyrs),
               seq(6, 6, len=nyrs),
               seq(7, 7, len=nyrs))

## Define the results lists to fill in loop below. There's two
## ways to test  these, first their bias for the total biomass,
## and second to what they are actually trying to estimate. The
## second is more of a self check.
index.combined.total.list <- index.combined.self.list <-
  index.ats.total.list <- index.ats.self.list <-
    index.bts.total.list <- index.bts.self.list <-
      pars.ats.list <- pars.bts.list <- pars.combined.list <- list()
kk <- 1
## Start of looping. Out loop is over trend in the OM, inner loop
## is replicates of the OM with process error.
for(trend in c('trend','flat')){
  if(trend=='trend'){
    beta1 <- beta1.trend; beta2 <- beta2.trend
  } else {
    beta1 <- beta1.flat; beta2 <- beta2.flat
  }
  ## Step 2: For each replicate build the "truth" and generate
  ## data from it
  par.truth <- par.truth0
  par.truth[grep('beta1_ft', par.names)] <-  as.vector(t(beta1))
  par.truth[grep('beta2_ft', par.names)] <-  as.vector(t(beta2))
  ## Loosely based on fitted model
  par.truth[grep('gamma1_ctp', par.names)] <- c(.4,.2,.5)
  par.truth[grep('gamma2_ctp', par.names)] <- c(.9, -.3, -.3)
  par.truth[grep('lambda1_k', par.names)] <- 0
  ## Half of the estimates from pollock just so we need fewer
  ## runs to see the patterns clearly
  par.truth[grep('L_omega1_z', par.names)] <- c(.7, 1.8, 1.9)/2
  par.truth[grep('L_omega2_z', par.names)] <- c(2, .5, .5)/2
  par.truth[grep('logSigmaM', par.names)] <- c(.4, .8)/2
  ## From pcod because we didn't estimate them for pollock
  par.truth[grep('ln_H_input', par.names)] <- c(.3235, -1.209)
  par.truth[grep('logkappa1', par.names)] <- log(sqrt(8)/400)
  par.truth[grep('logkappa2', par.names)] <- log(sqrt(8)/350)
  ## all random effecst are zero since resimulated internally anyway
  par.truth[Obj.OM$env$random] <- 0
  stopifnot(all(!is.na(par.truth)))
  ## tmp <- t(Obj.OM$report(par.truth)$Index_cyl[,,1])
  ## matplot(log(tmp), type='b')
  ## tmp <- t(apply(tmp, 1, function(x) cumsum(x)/sum(x)))
  ## matplot(tmp, type='b', ylim=c(0,1))
  for(iii in 1:nsim){
    Data_Geostat <- dat0
    set.seed(iii) # works with TMB?? probably not
    ## These are the truths after simulating new random effects
    ## (process error). Using the built-in TMB simulate feature.
    simdat <- Obj.OM$simulate(par=par.truth, complete=TRUE)
    index.total.truth <- log(simdat$ColeIndex_cy[1,])
    index.bts.truth <- log(simdat$ColeIndex_cy[2,])
    index.ats.truth <- log(simdat$ColeIndex_cy[3,])
    index.stratum1.truth <- log(simdat$Index_cyl[1,,1])
    index.stratum2.truth <- log(simdat$Index_cyl[2,,1])
    index.stratum3.truth <- log(simdat$Index_cyl[3,,1])
    ## Rebuild the Obj with the new simulated data
    Data_Geostat$Catch_KG <- simdat$b_i
    ## Process the simulated data to be put into the model
    DF1 <- subset(Data_Geostat, Gear=='BT')
    DF2 <- subset(Data_Geostat, Gear=='AT2')
    DF3 <- subset(Data_Geostat, Gear=='AT3')
    DF1$knot_i <- DF2$knot_i <- DF3$knot_i <- NULL
    ## Check that there are no 100% encounters or non-encounters
    ## for the combined
    pct.zero <- Data_Geostat %>% group_by(Year, Gear) %>%
      summarize(pct.zero=mean(Catch_KG==0)) %>% pull(pct.zero)
    ## AT is a special case since gets added together so recreate
    ## that and catch it here.
    pct.zero2 <- DF2 %>% mutate(Gear=factor('AT'), Catch=Catch_KG+DF3$Catch_KG) %>%
      group_by(Year) %>% summarize(pct.zero=mean(Catch==0)) %>% pull(pct.zero)
    if(any(pct.zero==0) | any(pct.zero2==0)) {
      warning(paste("100% encounter in", iii, "...skipping"))
      next
    }
    if(any(pct.zero==1) | any(pct.zero2==1)) {
      warning(paste("0% encounter in", iii, "...skipping"))
      next
    }
    ## First run the combined model
    control <- list(beta2temporal=TRUE, n_x=n_x, model='combined',
                    n_eps1=0, n_eps2=0, n_omega2="IID", n_omega1="IID",
                    beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                    temporal=0, fixlambda=2, aniso=TRUE, kappaoff=0,
                    simdata=TRUE, simulation=TRUE, make_plots=iii==1)
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_combined")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper,
                   ## start from MLE to speed things up
                   startpar=par.truth[-Obj.OM$env$random],
                   lower=TmbList$Lower, control=list(trace=0),
                   loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    clean.dir(savedir)
    if(iii==1) plot.vastfit(results, plotmaps=TRUE)
    index.combined.self.list[[kk]] <-
      data.frame(rep=iii, trend=trend, results$Index.strata,
                 maxgrad=Opt$max_gradient,
                 truth=as.vector(log(t(simdat$Index_cyl[,,1]))))
    index.combined.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=as.vector(log(t(simdat$ColeIndex_cy))))
    pars.combined.list[[kk]] <-
      data.frame(rep=iii, trend=trend, par=names(Opt$par),
                 par.num=1:length(Opt$par), maxgrad=Opt$max_gradient,
                 est=Opt$par, truth=par.truth[-Obj.OM$env$random])
    ## Now repeat with AT
    control$model <- 'ats'; control$make_plots <- FALSE
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_ats")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=0), loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    if(iii==1) plot.vastfit(results, plotmaps=TRUE)
    clean.dir(savedir)
    index.ats.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.total.truth)
    index.ats.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.ats.truth)
    ## Now repeat with BT
    control$model <- 'bts'
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_bts")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=0), loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    if(iii==1) plot.vastfit(results, plotmaps=TRUE)
    clean.dir(savedir)
    index.bts.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.total.truth)
    index.bts.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.bts.truth)
    kk <- kk+1
    ## Unlink DLLs to prevent error b/c there's a max # that can be loaded
  }
}

### Process results and save them to file
index.combined.self <- do.call(rbind, index.combined.self.list)
index.combined.total <- do.call(rbind, index.combined.total.list)
index.bts.self <- do.call(rbind, index.bts.self.list)
index.bts.total <- do.call(rbind, index.bts.total.list)
index.ats.self <- do.call(rbind, index.ats.self.list)
index.ats.total <- do.call(rbind, index.ats.total.list)
index.total <- rbind(index.ats.total, index.bts.total, filter(index.combined.total, strata=='total'))
index.self <-  rbind(index.ats.self, index.bts.self, filter(index.combined.self))
pars.combined <- do.call(rbind, pars.combined.list)
saveRDS(list(index.combined.self=index.combined.self,
             index.combined.total=index.combined.total,
             index.total=index.total, index.self=index.self,
             pars.combined=pars.combined),
        file='results/simulation.RDS')

x <- readRDS('results/simulation.RDS')

### Quick plots of these results
## Performance relative to total
g <- ggplot(x$index.total, aes(year, (est-truth)/truth, group=rep, color=log(maxgrad))) + geom_line() +
  facet_grid(trend~strata) + theme_bw() +
  geom_hline(yintercept=0, col=2)
ggsave('plots/simulation_RE_total.png', g, width=7, height=5)
g <- filter(x$index.self, model!='combined') %>%
  ggplot(aes(year, (est-truth)/truth, group=interaction(rep,strata), color=log(maxgrad))) +
  geom_line() +  geom_hline(yintercept=0, col=2) +
  geom_hline(yintercept=0, col=2) +
  facet_grid(trend~model) + theme_bw()
ggsave('plots/simulation_RE_self.png', g, width=7, height=5)
g <- filter(x$index.self, model=='combined') %>%
  ggplot(aes(year, (est-truth)/truth, group=interaction(rep,strata), color=log(maxgrad))) +
  geom_line() +  geom_hline(yintercept=0, col=2) +
  facet_grid(trend~strata) + theme_bw()
ggsave('plots/simulation_RE_self_combined.png', g, width=7, height=5)
g <- filter(x$pars.combined, !grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth)/truth)) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x') + theme_bw() +  geom_hline(yintercept=0, col=2)
ggsave('plots/simulation_RE_notbetas.png', g, width=7, height=5)
g <- filter(x$pars.combined, grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth))) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x') + theme_bw() +
  geom_hline(yintercept=0, col=2)
ggsave('plots/simulation_RE_betas.png', g, width=7, height=5)
out <- filter(x$pars.combined, grepl('beta', par)) %>%
  cbind(stratum=c(1,2,3), year=rep(1:8, each=3)) %>%
  mutate(stratum=factor(stratum, levels=1:3, labels=c('<0.5m', '0.5-16m', '>16')),
         abs.error=(est-truth)) %>%
  select(-par.num, -maxgrad, -est, -truth) %>%
  spread(key=par, value=abs.error)
g <- ggplot(out, aes(beta1_ft, beta2_ft, color=trend)) + geom_point() +
  facet_grid(year~stratum) + geom_hline(yintercept=0, color=2) +
  geom_vline(xintercept=0, color=2)
ggsave('plots/simulation_RE_betas_pairwise.png', g, width=12, height=9)
## Get a proportion of biomass across strata in each replicate
g <- x$index.combined.self %>%
  select(rep, trend, year, strata, truth) %>%
  group_by(rep, trend, year) %>%
  mutate(pct=(exp(truth))/sum(exp(truth))) %>%
  group_by(trend, year, strata) %>%
  summarize(mean.pct=mean(pct)) %>%
  ## spread(strata, mean.pct) %>%
  mutate(stratum=factor(strata, levels=c('stratum3', 'stratum2', 'stratum1'),
                        labels=rev(c('<0.5m', '0.5-16m',
                                     '>16m')))) %>%
  ggplot(aes(year, mean.pct, fill=stratum)) + geom_area() +
  facet_wrap('trend', nrow=2) + ylab('Proportion Abundance') + theme_bw()
ggsave('plots/simulation_OM_proportions.png', g, width=7, height=5)
## Look at the simulated truth
g <- ggplot(filter(x$index.self, model=='combined'), aes(year, y=truth, group=rep)) +
  geom_line() + facet_grid(trend~strata) + theme_bw()
ggsave('plots/simulation_OM.png', g, width=7, height=5)
g <- ggplot(filter(x$index.total, strata=='total'), aes(year, y=truth, group=rep)) +
  geom_line() + facet_grid(trend~strata) + theme_bw()
ggsave('plots/simulation_OM_total.png', g, width=7, height=5)




