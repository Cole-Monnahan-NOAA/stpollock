### This files runs the simulation testing component of the analysis.

nsim <- 50
n_x <- 100 # number of knots for OM and EM (they match)
ns <- 1 # number of newton steps
ln <- 5 # loop number in optimizer
clean.dir <- function(savedir){
  ## Unlink and delete TMB objects to save space and prevent errors about
  ## too many dlls loaded.
  dyn.unload(paste0(savedir,'/VAST_v8_0_0.dll'))
  x <- paste0(savedir,'/', c('VAST_v8_0_0.dll', 'VAST_v8_0_0.o', 'Record.RData'))
  trash <- file.remove(x)
}

### Step 1: Build a base OM model from which to simulate
### data. This model has no ST effects and independent spatial
### effects among strata. We use the anisotropy and kappa pars
### from the base case model and estimate them here b/c we use
### MLE.
control <- list(beta2temporal=TRUE, n_x=n_x, model='combined',
                n_eps1='IID', n_eps2='IID', n_omega2="IID", n_omega1='IID',
                beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                kappaoff=0, temporal=0, fixlambda=1,
                simdata=FALSE, simulation=TRUE,
                aniso=TRUE)
controlOM <- control
savedir <- paste0(getwd(), '/simulations/OM')
source("prepare_inputs.R")
Obj.OM <- Obj
par.names <- names(Obj.OM$env$last.par)
par.truth0 <- Obj.OM$env$last.par*NA ## make sure everything is changed
dat0 <- Data_Geostat ## save a copy of original data since gets overwritten
nyrs <- length(unique(Data_Geostat$Year))
beta1.trend <- cbind(seq(1,-1.5, len=nyrs),
               seq(.5,.5, len=nyrs),#c(seq(-1, .5, len=nyrs-3), .5,.5,.5),
               seq(-.5,1, len=nyrs))
beta2.trend <- cbind(seq(6,4, len=nyrs),
               seq(5, 5, len=nyrs),
               seq(4, 5.5, len=nyrs))
beta1.flat <- cbind(seq(1,1, len=nyrs),
               seq(-1, -1, len=nyrs),
               seq(-1,-1, len=nyrs))
beta2.flat <- cbind(seq(5,5, len=nyrs),
               seq(6, 6, len=nyrs),
               seq(7, 7, len=nyrs))

## Define the results lists to fill in loop below. There's two
## ways to test  these, first their bias for the total biomass,
## and second to what they are actually trying to estimate. The
## second is more of a self check.
indexc.total.list <- indexc.self.list <-
  indexa.total.list <- indexa.self.list <-
    indexb.total.list <- indexb.self.list <-
      pars.list <- list()
kk <- 1
## Start of looping. Outer loop is over trend in the OM, inner loop
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
  par.truth[grep('lambda2_k', par.names)] <- 0
  ## Halve and quarter the estimates from pollock just so we need
  ## fewer runs to see the patterns clearly
  par.truth[grep('L_omega1_z', par.names)] <- c(.7, 1.8, 1.9)/2
  par.truth[grep('L_omega2_z', par.names)] <- c(2, .5, .5)/2
  par.truth[grep('L_epsilon1_z', par.names)] <- c(.6,.85, 1.3)/4
  par.truth[grep('L_epsilon2_z', par.names)] <- c(1.3, 1.2, .9)/4
  par.truth[grep('logSigmaM', par.names)] <- 1000*c(.4, .8)/2
  ## From base case model
  par.truth[grep('ln_H_input', par.names)] <- c(.29, -.73)
  par.truth[grep('logkappa1', par.names)] <- -5.1
  par.truth[grep('logkappa2', par.names)] <- -4.9
  ## all random effecst are zero since resimulated internally anyway
  par.truth[Obj.OM$env$random] <- 0
  stopifnot(all(!is.na(par.truth)))
  ##  ## Make quick plot of OM without any process error
  ## png(paste0('plots/simulation_OM_', trend, '.png'), width=7,
  ##     height=5, units='in', res=500)
  ## par(mfrow=c(1,3), mar=c(4,4,1,1), mgp=c(2,.5,0), tck=-.02,
  ##     oma=c(0,0,3,0), cex.axis=.8, col.axis=gray(.3))
  ## tmp <- t(Obj.OM$report(par.truth)$Index_cyl[,,1])
  ## matplot(log(tmp), type='b', ylab='log-index')
  ## tmp2 <- cbind(tmp[,1]+tmp[,2], tmp[,2]+tmp[,3], rowSums(tmp))
  ## matplot(log(tmp2), type='b', ylab='log-index')
  ## tmp <- t(apply(tmp, 1, function(x) cumsum(x)/sum(x)))
  ## matplot(tmp, type='b', ylim=c(0,1), ylab='Percent density')
  ## mtext(paste('Simulation model for:', trend), outer=TRUE,
  ##       line=0, cex=1.5)
  ## mtext('Year', side=1, line=-2, outer=TRUE)
  ## dev.off()
  for(iii in 2:nsim){
    Data_Geostat <- dat0
    set.seed(iii) # works with TMB?? probably not
    ## These are the truths after simulating new random effects
    ## (process error). Using the built-in TMB simulate feature.
    simdat <- Obj.OM$simulate(par=par.truth, complete=TRUE)
    index.total.truth <- log(simdat$ColeIndex_cy[1,])
    indexb.truth <- log(simdat$ColeIndex_cy[2,])
    indexa.truth <- log(simdat$ColeIndex_cy[3,])
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
    control <- controlOM
    control$simdata <- TRUE; control$make_plots <- iii==1
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_combined/")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper,
                   ## start from MLE to speed things up
                   startpar=par.truth[-Obj.OM$env$random],
                   lower=TmbList$Lower, control=list(trace=5),
                   newtonsteps=ns,
                   loopnum=ln, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    clean.dir(savedir)
    if(iii==1)
      plot.vastfit(results, savedir, plotmaps=TRUE)
    indexc.self.list[[kk]] <-
      data.frame(rep=iii, trend=trend, results$Index.strata,
                 maxgrad=Opt$max_gradient,
                 truth=as.vector(log(t(simdat$Index_cyl[,,1]))))
    indexc.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=as.vector(log(t(simdat$ColeIndex_cy))))
    pars.list[[kk]] <-
      data.frame(rep=iii, trend=trend, par=names(Opt$par),
                 par.num=1:length(Opt$par), maxgrad=Opt$max_gradient,
                 est=Opt$par, truth=par.truth[-Obj.OM$env$random])
    ## Now repeat with AT
    control$model <- 'ats'; control$make_plots <- FALSE
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_ats/")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=5), loopnum=ln,
                   newtonsteps=ns, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    if(iii==1) plot.vastfit(results, savedir, plotmaps=TRUE)
    clean.dir(savedir)
    indexa.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.total.truth)
    indexa.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=indexa.truth)
    ## Now repeat with BT
    control$model <- 'bts'
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_bts/")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=5), loopnum=ln,
                   newtonsteps=ns, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    if(iii==1) plot.vastfit(results, savedir, plotmaps=TRUE)
    indexb.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=index.total.truth)
    indexb.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient, truth=indexb.truth)
    kk <- kk+1
    ## Unlink DLLs to prevent error b/c there's a max # that can be loaded
    clean.dir(savedir)
    ## Save temporary file in case it crashes
    save.image('simulations/.RData') #
  }
}

### Process results and save them to file
indexc.self <- do.call(rbind, indexc.self.list)
indexc.total <- do.call(rbind, indexc.total.list)
indexb.self <- do.call(rbind, indexb.self.list)
indexb.total <- do.call(rbind, indexb.total.list)
indexa.self <- do.call(rbind, indexa.self.list)
indexa.total <- do.call(rbind, indexa.total.list)
index.total <- rbind(indexa.total, indexb.total, dplyr::filter(indexc.total, strata=='total'))
index.self <-  rbind(indexa.self, indexb.self, indexc.self)
index.self$strata <- factor(index.self$strata,
                            levels=c('ats', 'bts', 'stratum1', 'stratum2', 'stratum3'),
                            labels=c('AT', 'BT', '<0.5m', '0.5-16m', '>16m'))
index.total$strata <- factor(index.total$strata,
                            levels=c('ats', 'bts', 'total'),
                            labels=c('AT', 'BT', 'Total'))
index.total$model <- factor(index.total$model,
                            levels=c('ats', 'bts', 'combined'),
                            labels=c('AT', 'BT', 'Combined'))
index.self$model <- factor(index.self$model,
                            levels=c('ats', 'bts', 'combined'),
                            labels=c('AT', 'BT', 'Combined'))
pars <- do.call(rbind, pars.list)
saveRDS(list(indexc.self=indexc.self, indexc.total=indexc.total,
             index.total=index.total, index.self=index.self,
             pars=pars), file='results/simulation.RDS')

x <- readRDS('results/simulation.RDS')
meta <- filter(x$index.self, year ==1 & !strata %in% c('<0.5m', '>16m'))
table.simulation <-  meta %>% group_by(model, trend) %>%
  dplyr::summarize(n=n(), pct.badgrads=mean(maxgrad>.001))
write.csv('results/table.simulation.csv', x=table.simulation)

### Quick plots of these results
library(cowplot)
theme_set(theme_bw())
alpha <- .5
mylim <- coord_cartesian(ylim=c(-.25,.25))
g <- ggplot(meta, aes(x=trend, y=maxgrad)) +
  geom_violin() + facet_wrap('model') +
  geom_hline(yintercept=(.001), col=2) + scale_y_log10()
ggsave('plots/simulation_maxgrads.png', g, width=7, height=5)
## Filter out the unconverged ones
indexc.self <- filter(x$indexc.self, maxgrad<=.001)
indexc.total <- filter(x$indexc.total, maxgrad<=.001)
index.self <- filter(x$index.self, maxgrad<=.001)
index.total <- filter(x$index.total, maxgrad<=.001)
pars <- filter(x$pars, maxgrad<=.001)

## Performance relative to total
g1 <- filter(index.self, model=='Combined') %>%
  ggplot(aes(factor(year), (est-truth)/truth)) +
  geom_violin() +  geom_hline(yintercept=0, col=2) +
  facet_grid(trend~strata) + theme_bw() + mylim +
  ylab('Relative Error')  + xlab('year')
g2 <- filter(index.self, model!='Combined') %>%
  ggplot(aes(factor(year), (est-truth)/truth)) +
  geom_violin() +  geom_hline(yintercept=0, col=2) +
  geom_hline(yintercept=0, col=2) +
  facet_grid(trend~model) + theme_bw() + mylim +
  ylab('Relative Error') + xlab('year')
g <- plot_grid(g1, g2, labels = c('A', 'B'), label_size = 12, nrow=2)
ggsave('plots/simulation_self_RE.png', g, width=7, height=9)

g <- ggplot(index.total, aes(factor(year), (est-truth)/truth)) +
 geom_violin(fill=gray(.9), scale='width') +
  geom_boxplot(width=.2, color=gray(.3), outlier.color=1, outlier.size=1) +
  facet_grid(trend~model) +
  theme_bw() +
  geom_hline(yintercept=0, col=2) +
  ylab('Error relative to total biomass') + xlab('Year')
ggsave('plots/simulation_total_RE.png', g, width=7, height=5)

g1 <- filter(pars, !grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth)/truth)) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x') + theme_bw() +  geom_hline(yintercept=0, col=2)
g2 <- filter(pars, grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth))) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x') + theme_bw() +
  geom_hline(yintercept=0, col=2)
g <- plot_grid(g1, g2, labels = c('A', 'B'), label_size = 12, nrow=2)
ggsave('plots/simulation_pars_RE.png', g, width=7, height=9)
## out <- filter(pars, grepl('beta', par)) %>%
##   cbind(stratum=c(1,2,3), year=rep(1:8, each=3)) %>%
##   mutate(stratum=factor(stratum, levels=1:3, labels=c('<0.5m', '0.5-16m', '>16')),
##          abs.error=(est-truth)) %>%
##   select(-par.num, -maxgrad, -est, -truth) %>%
##   spread(key=par, value=abs.error)
## g <- ggplot(out, aes(beta1_ft, beta2_ft, color=trend)) + geom_point() +
##   facet_grid(year~stratum) + geom_hline(yintercept=0, color=2) +
##   geom_vline(xintercept=0, color=2)
## ggsave('plots/simulation_RE_betas_pairwise.png', g, width=12, height=9)
## Get a proportion of biomass across strata in each replicate
g <- indexc.self %>%
  select(rep, trend, year, strata, truth) %>%
  group_by(rep, trend, year) %>%
  mutate(pct=(exp(truth))/sum(exp(truth))) %>%
  group_by(trend, year, strata) %>%
  summarize(mean.pct=mean(pct)) %>%
  ## spread(strata, mean.pct) %>%
  dplyr::mutate(stratum=factor(strata, levels=c('stratum3', 'stratum2', 'stratum1'),
                        labels=rev(c('<0.5m', '0.5-16m', '>16m')))) %>%
  ggplot(aes(year, mean.pct, fill=stratum)) + geom_area() +
  facet_wrap('trend', nrow=2) + ylab('Proportion Abundance')
ggsave('plots/simulation_OM_proportions.png', g, width=7, height=5)

## Look at the simulated truth by strata
g1 <- filter(index.self, model=='Combined') %>%
  ggplot(aes(year, y=truth, group=rep)) +
  geom_line(alpha=alpha) + facet_grid(trend~strata) +
  theme_bw() + labs(x='Year', y='log-index (truth)')
## Truth by gear
g2 <- rbind(filter(index.total, strata=='Total'),
           cbind(filter(index.self, model!='Combined'))) %>%
  ggplot(aes(year, truth, group=rep)) + geom_line(alpha=alpha) +
  facet_grid(trend~strata) + theme_bw() +
 labs(x='Year', y='log-index (truth)')
g <- plot_grid(g1, g2, labels = c('(a)', '(b)'), label_size = 12, nrow=2)
ggsave('plots/simulation_OM_truth.png', g, width=7, height=7)




