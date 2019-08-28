### This files runs the simulation testing component of the analysis.

rm(list=ls())
source("startup.R")
chains <- 4
options(mc.cores = chains)
nsim <- 50

## Build a base OM model from which to simulate data.
control <- list(beta2temporal=TRUE, n_x=100, model='combined',
                n_eps1=0, n_eps2=0, n_omega2=0, n_omega1=0,
                beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                kappaoff=12, temporal=0, fixlambda=2)
control$simdata <- FALSE
control$simulation <- TRUE
savedir <- paste0(getwd(), '/simulations/simfit_OM')
source("prepare_inputs.R")
Obj.OM <- Obj
par.names <- names(par)
par.truth0 <- par*NA ## make sure everything is changed
dat0 <- Data_Geostat ## save a copy of original data since gets overwritten
nyrs <- length(unique(Data_Geostat$Year))
beta1.trend <- cbind(seq(.2,-2, len=nyrs),
               seq(-2, -1, len=nyrs),
               seq(-.5,.1, len=nyrs))
beta2.trend <- cbind(seq(2,2.5, len=nyrs),
               seq(.5, 2, len=nyrs),
               seq(0, 1, len=nyrs))
beta1.flat <- cbind(seq(.2,.2, len=nyrs),
               seq(-1, -1, len=nyrs),
               seq(.1,.1, len=nyrs))
beta2.flat <- cbind(seq(2,2, len=nyrs),
               seq(.5, .5, len=nyrs),
               seq(1, 1, len=nyrs))


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
for(trend in rev(c('flat', 'trend'))){
  if(trend=='trend'){
    beta1 <- beta1.trend; beta2 <- beta2.trend
  } else {
    beta1 <- beta1.flat; beta2 <- beta2.flat
  }
  ## Build the OM true parameters.
  par.truth <- par.truth0
  par.truth[grep('beta1_ft', par.names)] <-  as.vector(t(beta1))
  par.truth[grep('beta2_ft', par.names)] <-  as.vector(t(beta2))
  par.truth[grep('gamma1_ctp', par.names)] <- -.1
  par.truth[grep('gamma2_ctp', par.names)] <- -.1
  par.truth[grep('lambda1_k', par.names)] <- -.05
  par.truth[grep('logSigmaM', par.names)] <- c(1, 1)
  ## matplot(t(rep.truth$Index_cyl[,,1]))
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
    ## Check that there are no 100% encounters or non-encounters
    pct.zero <- Data_Geostat %>% group_by(Year, Gear) %>%
      summarize(pct.zero=mean(Catch_KG==0)) %>% pull(pct.zero)
    if(any(pct.zero==0)) {
      warning(paste("100% encounter in", iii, "...skipping"))
      break
    }
    if(any(pct.zero==1)) {
      warning(paste("0% encounter in", iii, "...skipping"))
      break
    }
    ## Process the simulated data to be put into the model
    DF1 <- subset(Data_Geostat, Gear=='BT')
    DF2 <- subset(Data_Geostat, Gear=='AT2')
    DF3 <- subset(Data_Geostat, Gear=='AT3')
    DF1$knot_i <- DF2$knot_i <- DF3$knot_i <- NULL
    ## First run the combined model
    control <- list(beta2temporal=TRUE, n_x=100, model='combined',
                    n_eps1=0, n_eps2=0, n_omega2=0, n_omega1=0,
                    beta1temporal=TRUE, filteryears=TRUE, finescale=FALSE,
                    kappaoff=12, temporal=0, fixlambda=2,
                    simdata=TRUE, simulation=TRUE)
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_simfit_combined")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, startpar=par.truth,
                   lower=TmbList$Lower, control=list(trace=10),
                   newtonsteps=1, loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    ## fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
    ##                iter=800, open_progress=FALSE, warmup=700,
    ##                init='last.par.best', thin=1,
    ##                control=list(max_treedepth=5))
    ## plot.vastfit(results, plotmaps=TRUE)
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
                 par.num=1:length(Opt$par),
                 maxgrad=Opt$max_gradient,
                 est=Opt$par, truth=par.truth)
    ## Now repeat with AT
    control$model <- 'ats'
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_simfit_bts")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=0), newtonsteps=1,
                   loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    ## fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
    ##                iter=800, open_progress=FALSE, warmup=700,
    ##                init='last.par.best', thin=1,
    ##                control=list(max_treedepth=5))
    ## plot.vastfit(results, plotmaps=TRUE)
    index.ats.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=index.total.truth)
    index.ats.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=index.ats.truth)
    ## Now repeat with BT
    control$model <- 'bts'
    savedir <- paste0(getwd(), '/simulations/', trend, "_", iii, "_simfit_ats")
    source("prepare_inputs.R")
    Opt <- fit_tmb(TmbList$Obj, upper=TmbList$Upper, lower=TmbList$Lower,
                   control=list(trace=0), newtonsteps=1,
                   loopnum=5, getsd=FALSE)
    results <- process.results(Opt, Obj, Inputs, model, space, savedir)
    ## fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
    ##                iter=800, open_progress=FALSE, warmup=700,
    ##                init='last.par.best', thin=1,
    ##                control=list(max_treedepth=5))
    ## plot.vastfit(results, plotmaps=TRUE)
    index.bts.total.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=index.total.truth)
    index.bts.self.list[[kk]] <-
      data.frame(rep=iii,results$Index, trend=trend,
                 maxgrad=Opt$max_gradient,
                 truth=index.bts.truth)
    kk <- kk+1
  }
}
index.combined.self <- do.call(rbind, index.combined.self.list)
index.combined.total <- do.call(rbind, index.combined.total.list)
index.bts.self <- do.call(rbind, index.bts.self.list)
index.bts.total <- do.call(rbind, index.bts.total.list)
index.ats.self <- do.call(rbind, index.ats.self.list)
index.ats.total <- do.call(rbind, index.ats.total.list)
index.total <- rbind(index.ats.total, index.bts.total, filter(index.combined.total, strata=='total'))
index.self <-  rbind(index.ats.self, index.bts.self, filter(index.combined.self))
pars.combined <- do.call(rbind, pars.combined.list)


## Look at the simulated truth
ggplot(filter(index.self, model=='combined'), aes(year, y=truth, group=rep)) +
  geom_line() + facet_grid(trend~strata)
## Performance relative to total
ggplot(index.total, aes(year, (est-truth)/truth, group=rep, color=log(maxgrad))) + geom_line() +
  facet_grid(trend~strata)

filter(index.self, model!='combined') %>%
  ggplot(aes(year, (est-truth)/truth, group=interaction(rep,strata), color=log(maxgrad))) +
  geom_line() +  geom_hline(yintercept=0, col=2) + facet_grid(trend~model)
filter(index.self, model=='combined') %>%
  ggplot(aes(year, (est-truth)/truth, group=interaction(rep,strata), color=log(maxgrad))) +
  geom_line() +  geom_hline(yintercept=0, col=2) + facet_grid(trend~strata)

filter(pars.combined, !grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth))) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x')

filter(pars.combined, grepl('beta', par)) %>%
  ggplot(aes(factor(par.num), (est-truth)/truth)) +
  geom_violin() + geom_abline(slope=0, intercept=0, color='red')+
  facet_grid(trend~par, scales='free_x')

