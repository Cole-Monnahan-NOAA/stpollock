
library(TMB)
## library(devtools)
## remove.packages('VAST')
## install_github('James-Thorson/FishStatsUtils', ref='2.3.2')
## install_github('James-Thorson/VAST', ref='3.2.1')

library(reshape2)
## library(plyr)
library(TMBhelper)
library(snowfall)
library(abind)
library(tmbstan)
library(shinystan)
## library(magrittr)
library(tidyverse)
## library(here)
library(VAST)
library(FishStatsUtils)
library(maps)
library(mapdata)
library(cowplot)
Version <- "VAST_v8_0_0"
compile('models/VAST_v8_0_0.cpp')

strata.labels.combined <- c('<0.5m', '0.5-16m', '>16m')
strata.labels.ats <- '>0.5m'
strata.labels.bts <- '<16m'


dir.create('simulations', showWarnings=FALSE)
dir.create('plots', showWarnings=FALSE)
dir.create('simulations/plots', showWarnings=FALSE)

## source("simulator.R")

## This generic function uses the Obj on the global space to generate
## random inits by jittering the defaults.
prior.fn <- function(){
  par.all <- Obj$env$last.par
  fixed <- par.all[-Obj$env$random]
  random <- par.all[Obj$env$random]
  fixed[fixed==0] <- 1
  fixed <- fixed*runif(length(fixed), min=.5, max=1.5)
  random <- rnorm(length(random),0,.1)
  par.all[-Obj$env$random] <- fixed
  par.all[Obj$env$random] <- 0*random
  par.all
}


plot.change <- function(Report){
  ## Look at trend in % of population <3m.
  Dtmp <- Report$D_gcy
  dimnames(Dtmp) <- list(cell=1:control$n_x, stratum=strata.labels.combined, year=years)
  Dtmp.wide <- plyr::melt(Dtmp)
  ## pct <- ddply(Dtmp.wide, .(cell, year), mutate, pct=value/sum(value))
  pct <- Dtmp.wide %>% group_by(cell, year) %>% mutate(pct=value/sum(value))
  ## g <- ggplot(pct, aes(year, pct, fill=stratum)) + geom_area() + facet_wrap('cell')
  pct2 <- subset(pct, stratum==strata.labels.combined[1])
  MatDat <- plyr::ddply(pct2, .(cell), function(x) {
    fit <- lm(pct~year, data=x)
    out <- data.frame(slope=coef(fit)[2], pvalue=summary(fit)$coefficients[2,4] )
    row.names(out) <- NULL
    out
  })
  temp <- merge(pct2, MatDat, by='cell')
  temp$significant <- temp$pvalue<.05
  ## g <- ggplot(temp, aes(year, pct, group=cell, color=slope)) + geom_line() + facet_wrap('significant')
  ## MatDat[,2] <- ifelse(MatDat$pvalue<.05, MatDat[,2], 0)
  mdl <- make_map_info(Region=Region, spatial_list=Spatial_List,
                       Extrapolation_List=Extrapolation_List )
  PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(MatDat[,2], ncol=1),
             PlotDF=mdl$PlotDF, zlim=range(MatDat[,2]),
             MapSizeRatio=mdl$MapSizeRatio,
             Xlim=mdl$Xlim, Ylim=mdl$Ylim,
             FileName=paste0(savedir, '/map_change'),
             Year_Set=Year_Set[1], Legend=mdl$Legend,
             mfrow = c(1,1),
             zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
             oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)

  ## Look at correlation among strata in each grid cell
  Dtmp <- Report$Epsilon1_gct
  if(!is.null(Dtmp)){
    dimnames(Dtmp) <- list(cell=1:control$n_x, stratum=strata.labels.combined, year=years)
    Dtmp.wide <- plyr::melt(Dtmp)
    Dtmp.long <- dcast(Dtmp.wide, year+cell~stratum, value.var='value')
    names(Dtmp.long)[3:5] <- strata.labels.combined
    g <- ggplot(Dtmp.long, aes(x=stratum1, y=stratum2, color=year)) +
      geom_point()+ facet_wrap('cell') +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ggsave(file.path(savedir, 'eps1_correlations_pairwise.png'), g, width=12, height=9)
    cors <- plyr::ddply(Dtmp.wide, .(cell), function(x){
      s1=subset(x, stratum==strata.labels.combined[1])$value
      s2=subset(x, stratum==strata.labels.combined[2])$value
      s3=subset(x, stratum==strata.labels.combined[3])$value
      data.frame(cor12=cor(s1,s2), cor13=cor(s1,s3), cor23=cor(s2,s3))
    })
    cors.wide <- plyr::melt(cors[,-1], measure.vars=c('cor12', 'cor13', 'cor23'))
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,2], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps1_cor12'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,3], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps1_cor13'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,4], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps1_cor23'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
  }
  ## Repeat for eps2
  Dtmp <- Report$Epsilon2_gct
  if(!is.null(Dtmp)){
    dimnames(Dtmp) <- list(cell=1:control$n_x, stratum=strata.labels.combined, year=years)
    Dtmp.wide <- plyr::melt(Dtmp)
    Dtmp.long <- dcast(Dtmp.wide, year+cell~stratum, value.var='value')
    names(Dtmp.long)[3:5] <- strata.labels.combined
    g <- ggplot(Dtmp.long, aes(x=stratum1, y=stratum2, color=year)) +
      geom_point()+ facet_wrap('cell') +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ggsave(file.path(savedir, 'eps2_correlations_pairwise.png'), g, width=12, height=9)
    cors <- plyr::ddply(Dtmp.wide, .(cell), function(x){
      s1=subset(x, stratum==strata.labels.combined[1])$value
      s2=subset(x, stratum==strata.labels.combined[2])$value
      s3=subset(x, stratum==strata.labels.combined[3])$value
      data.frame(cor12=cor(s1,s2), cor13=cor(s1,s3), cor23=cor(s2,s3))
    })
    cors.wide2 <- plyr::melt(cors[,-1], measure.vars=c('cor12', 'cor13', 'cor23'))
    cors.wide3 <- rbind(cbind(variable='epsilon1',cors.wide), cbind(variable='epsilon2',cors.wide2))
    names(cors.wide3)[2] <- 'correlation'
    g <- ggplot(cors.wide3, aes(x=value)) + geom_histogram(position='identity', bins=20) +
      facet_grid(variable~correlation)
    ggsave(file.path(savedir, 'eps_correlations.png'), g, width=7, height=5)
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,2], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps2_cor12'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,3], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps2_cor13'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
    PlotMap_Fn(MappingDetails=mdl$MappingDetails, Mat=matrix(cors[,4], ncol=1),
               PlotDF=mdl$PlotDF, zlim=c(-1,1),
               MapSizeRatio=mdl$MapSizeRatio,
               Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_eps2_cor23'),
               Year_Set=Year_Set[1], Legend=mdl$Legend,
               mfrow = c(1,1),
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=.5, plot_legend_fig=FALSE, pch=16)
  }
}

get.resids.tmp <- function(case){
  type <- 'combined'
  if(length(grep('combinedoff', x=savedir))>0) type <- 'combinedoff'
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  sigtmp <- results$Report$SigmaM[as.numeric(Data_Geostat$Gear)]^2/2
  df <- data.frame(obs=log(Data_Geostat$Catch_KG),
                   predicted= log(results$Report$R2_i)-sigtmp,
                   gear=Data_Geostat$Gear,
                   year=Data_Geostat$Year, case=case, type=type)
  df <- subset(df, obs>0) ## drop zeroes
  df
}

get.index.tmp <- function(case){
  type <- 'combined'
  if(length(grep('combinedoff', x=savedir))>0) type <- 'combinedoff'
  results <- process.results(Opt, Obj, Inputs, model, space, savedir)
  index.model <- results$Index.strata
  if(type=='combined'){
    ## Need to recalculate the gears for BTS, ats1 and ats2
    tmp <- dcast(index.model[,c('year', 'strata', 'est')],
                 year~strata, value.var='est')
    tmp[,'BTS'] <- log(exp(tmp$stratum1)+exp(tmp$stratum2))
    tmp[,'ATS1'] <- tmp$stratum2
    tmp[,'ATS2'] <- tmp$stratum3
    index.model <- plyr::melt(tmp[, c('year', 'BT', 'AT2', 'AT3')],
                        id.vars='year', value.name='logdensity', variable.name='gear')
  } else {
    ## strata here are actually gears b/c not summing
    levels(index.model$strata) <- c("BT", 'AT2', 'AT3')
    index.model$gear <- index.model$strata
    index.model$logdensity <- index.model$est
  }
  index.model$model <- 'model'
  index.raw$model <- 'raw data'
  index.model$case <- index.raw$case <- case
  index.model$type <- index.raw$type <- type
  index.raw$logdensity <- log(index.raw$value)
  x <- c('case', 'model', 'type', 'year', 'gear', 'logdensity')
  out <- rbind(index.model[,x], index.raw[,x])
  out
}


locate.year <- function(yr){
  library(maps)
  library(geosphere)
  map('world2', xlim=c(180, 202), ylim=c(53,63), add=FALSE)
  map.axes()
  with(subset(ats, year==yr), points(360+lon,lat, pch=16, cex=.1))
  with(subset(bts, year==yr), points(360+lon, lat, pch=1, cex=.1, col=2))
  locs <- as.data.frame(locator())
  points(locs, pch=15)
  iseq <- seq(1, nrow(locs), by=2)
  out <- list()
  for(i in iseq){
    ## this gap is 1/2 nautical miles
    out[[i]] <- add.zeroes(locs[i,], locs[i+1,], gap=0.926)
    points(out[[i]], pch=5, cex=.1)
  }
  out <- do.call(rbind, out)
  df <- data.frame(X=-999, year=yr, lon=out[,1]-360, lat=out[,2],
                   dist=NA, surface=NA, ground=NA,
                   strata1=0, strata2=0, strata3=0)
  return(df)
}

## out <- add.zeroes(c(-175,59), c(-174, 56), gap=.01)
## plot(out, pch=5)
add.zeroes <- function(p1, p2, gap){
  ## Function to add points with zero observations between coordinates p1
  ## and p2 with a distance between them of gap (in km).
  p1 <- as.numeric(p1)
  p2 <- as.numeric(p2)
  ## Need to get distance on globe then calculate how many points to create
  ## equally space ones with the right gap.
  ## distance in km
  dist <- distHaversine(p1=c(p1[1]-360, p1[2]), p2=c(p2[1]-360, p2[2]))/1000
  steps <- round(dist/gap)
  if(p1[1]==p2[1]) p2[1] <- p2[1]+.01 # prevent infinite slope
  z2 <- p2-p1; z1 <- c(0,0)
  slope <- z2[2]/z2[1]
  if(slope<0)  xseq <- seq(z2[1], z1[1], length=steps)
  else xseq <- seq(z1[1], z2[1], length=steps)
  yseq <- xseq*slope
  out <- data.frame(lon=xseq+p1[1], lat=yseq+p1[2])
  out
}

process.results <- function(Opt, Obj, Inputs, model, space, savedir){
  Report  <-  Obj$report()
  ParHatList <- Obj$env$parList(Opt$par)
  ParHat <- Opt$par
  par2 <- do.call(c, sapply(unique(names(ParHat)), function(x){
    n <- length(which(names(ParHat)==x))
    if(n>1) paste0(x,"_",1:n) else x
  }))
  if(is.null(Opt$SD$cov.fixed)){ # no SD
    SE <- rep(Inf, len=length(ParHat))
  } else {
    SE <- sqrt(diag(Opt$SD$cov.fixed))
  }
  est <- data.frame(par=names(ParHat), par2=par2, est=ParHat, lwr=ParHat-1.96*SE,
                    upr=ParHat+1.96*SE, SE=SE)
  est$significant <- !(est$lwr<0 & est$upr>0)
  write.csv(est, file=paste0(savedir, "/estimates.csv"))
  Index <- calculate.index(Opt, Report, model, space, log=TRUE, strata=FALSE)
  Index.strata <- calculate.index(Opt, Report, model, space, log=TRUE, strata=TRUE)
  Save  <-  list(Index=Index, Opt=Opt, Report=Report, ParHat=ParHat,
                 ParHatList=ParHatList, est=est, Index.strata=Index.strata,
                 SE=SE, Inputs=Inputs, savedir=savedir, model=model)
  saveRDS(Save, file=paste0(savedir, '/Save.RDS'))
  return(Save)
}


plot.sampler.params <- function(fit){
  sp <- get_sampler_params(fit)
  sp <- lapply(1:length(sp), function(i) data.frame(chain=i, iter=1:nrow(sp[[i]]), sp[[i]])) %>% do.call(rbind,.) %>%
    as.data.frame() %>%
    mutate(log_stepsize=log(stepsize__), chain=factor(chain)) %>%
    select(-energy__, -n_leapfrog__, -stepsize__) %>%
    gather(variable, value, -chain, -iter) %>% filter(iter>5)
  g <- ggplot(sp, aes(iter, y=value, color=chain)) + geom_point(alpha=.5) +
    facet_wrap('variable', scales='free_y', ncol=1) + theme_bw()
  ggsave(paste0(savedir, '/sampler_params.png'), g, width=7, height=7)
}

get.results.mcmc <- function(Obj, fit){## Get parameters and drop log-posterior
  df <- as.matrix(fit)
  df <- df[,-ncol(df)] # drop lp__ column
  plot.sampler.params(fit)
  if(model=='combined'){
    strata <- strata.labels.combined
    gear <- c('Total', 'BT', 'AT')
  } else if(model=='ats'){
    strata <- strata.labels.ats
    gear <- 'AT'
  } else if(model=='bts'){
    strata <- strata.labels.bts
    gear <- 'BT'
  } else {
    stop("invalid model")
  }
  ## Save text file of results
  fixed <- (df[, -Obj$env$random])
  fixed.summary <- do.call(rbind, apply(fixed, 2, function(x)
    round(data.frame(median=median(x), mean=mean(x), sd=sd(x), min=min(x),
                     max=max(x), lwr=quantile(x, .025), upr=quantile(x,.975)),4)))
  fixed.summary <- cbind(par=rownames(fixed.summary), fixed.summary)
  write.csv(file=paste0(savedir,'/fixed.estimates.csv'), x=fixed.summary, row.names=FALSE)
  random <- (df[, Obj$env$random])
  random.summary <- do.call(rbind, apply(random, 2, function(x)
    round(data.frame(median=median(x), mean=mean(x), sd=sd(x), min=min(x),
                     max=max(x)),3)))
  random.summary <- cbind(par=rownames(random.summary), random.summary)
  write.csv(file=paste0(savedir,'/random.estimates.csv'), x=random.summary, row.names=FALSE)
  index.gear.tmp <- index.strata.tmp <- D_gcy.list <- list()
  covcor_omega1.list <- covcor_omega2.list <- covcor_epsilon1.list <- list()
  message("Looping through and calculating report...")
  tmp <- Obj$report(df[1,])
  ## Merge these into 4d arrays, last dimension is posterior draw number
  D_gcyn <- array(NA, dim=c(dim(tmp$D_gcy), nrow(df)))
  R1_in <- R2_in <- array(NA, dim=c(length(tmp$R1_i), nrow(df)))
  PR1_in <- PR2_in <- R1_in
  beta1_tcn <- beta2_tcn <- array(NA, dim=c(dim(tmp$beta1_tc), nrow(df)))
  R1_gcy.list <- R2_gcy.list <- list()
  for(i in 1:nrow(df)){
    if(i %% 50 ==0) print(i)
    tmp <- Obj$report(df[i,])
    beta1_tcn[,,i] <- tmp$beta1_tc
    beta2_tcn[,,i] <- tmp$beta2_tc
    index.strata.tmp[[i]] <-
      data.frame(year=rep(years, each=length(strata)), iter=i, density=log(as.numeric(tmp$Index_cy)),
                 stratum=strata)
    index.gear.tmp[[i]] <-
      data.frame(year=rep(years, each=length(gear)), iter=i, density=log(as.numeric(tmp$ColeIndex_cy)),
                 gear=gear)
    D_gcyn[,,,i] <- tmp$D_gcy
    R1_in[,i] <- tmp$R1_i
    R2_in[,i] <- tmp$R2_i
    R1_gcy.list[[i]] <- tmp$R1_gcy
    R2_gcy.list[[i]] <- tmp$R2_gcy
    covcor_omega1.list[[i]] <- tmp$lowercov_uppercor_omega1
    covcor_omega2.list[[i]] <- tmp$lowercov_uppercor_omega2
    covcor_epsilon1.list[[i]] <- tmp$lowercov_uppercor_epsilon1
    ## Calculate Pearson residuals
    for(ii in 1:nrow(PR1_in)){
      ## bernoulli for presence
      mui <- tmp$R1_i[ii]
      obs <- as.numeric(Data_Geostat$Catch_KG[ii]>0)
      PR1_in[ii,i] <- (obs-mui)/sqrt(mui*(1-mui)/1)
      ## log-normal for catch rate; NA for 0 observations
      obs <- Data_Geostat$Catch_KG[ii]
      if(obs>0){
        ## make sure to use the right variance as this depends on gear type
        gr <- as.numeric(Data_Geostat$Gear[ii])
        PR2_in[ii,i] <- (log(obs)-log(tmp$R2_i[ii])+tmp$SigmaM[gr]^2/2)/tmp$SigmaM[gr]
      }
    }
  }
  ## For now taking the mean of the Pearson residuals across posterior
  ## draws
  PR1 <- apply(PR1_in, 1, mean, na.rm=TRUE)
  PR2 <- apply(PR2_in, 1, mean, na.rm=TRUE)
  PR2[is.nan(PR2)] <- NA
  ## stopifnot(all.equal(D_gcyn[,,,1],D_gcy.list[[1]]))
  R1_gcyn <- array(do.call(c, R1_gcy.list), dim=c(dim(tmp$R1_gcy), nrow(df)))
  ## stopifnot(all.equal(R1_gcyn[,,,1],R1_gcy.list[[1]]))
  R2_gcyn <- array(do.call(c, R2_gcy.list), dim=c(dim(tmp$R2_gcy), nrow(df)))
  ## stopifnot(all.equal(R2_gcyn[,,,1],R2_gcy.list[[1]]))
  ## Report only the mean of these
  R1_gcy <- apply(R1_gcyn, 1:3, mean)
  R2_gcy <- apply(R2_gcyn, 1:3, mean)
  rm(R1_gcyn, R2_gcyn); gc()
  ## Organize the corcov matrices
  covcor_omega1 <- covcor_omega2 <- covcor_epsilon1 <- NULL
  if(length(covcor_omega1.list)>0)
    covcor_omega1 <- array(do.call(c, covcor_omega1.list), dim=c(3,3, nrow(df)))
  if(length(covcor_omega2.list)>0)
    covcor_omega2 <- array(do.call(c, covcor_omega2.list), dim=c(3,3, nrow(df)))
  if(length(covcor_epsilon1.list)>0)
    covcor_epsilon1 <- array(do.call(c, covcor_epsilon1.list), dim=c(3,3, nrow(df)))
  covcor <- list(covcor_omega1=covcor_omega1, covcor_omega2=covcor_omega2,
                 covcor_epsilon1=covcor_epsilon1)
  index.gear <- do.call(rbind, index.gear.tmp)
  index.strata <- do.call(rbind, index.strata.tmp)
  index.gear$gear <- as.factor(index.gear$gear)
  index.strata$stratum <- factor(index.strata$stratum, levels=strata)
  ## handle NaN's in density to prevent error and keep running other scenarios
  if(!all(is.finite(index.gear$density))){
    return(NULL)
  } else {
    index.gear2 <- plyr::ddply(index.gear, .(year, gear), summarize,
                         lwr=quantile(density, probs=.025),
                         upr=quantile(density, probs=.975),
                         est=median(density))
  }
  if(!all(is.finite(index.strata$density))){
    return(NULL)
  } else {
    index.strata2 <- plyr::ddply(index.strata, .(year, stratum), summarize,
                           lwr=quantile(density, probs=.025),
                           upr=quantile(density, probs=.975),
                           est=median(density))
  }
  ## Availability is in natural space and only makes sense for combined
  ## model
  if(model=='combined'){
    ## Massage to get the catchability by gear type
    tmp <- dcast(index.gear, year+iter~gear, value.var='density')
    availability <- within(tmp, {BT=exp(BT)/exp(Total);
      AT=exp(AT)/exp(Total)})
    availability <- plyr::melt(availability, id.vars=c('year', 'iter'),
                         measure.vars=c('AT', 'BT'),
                         variable.name='gear', value.name='availability')
    availability2 <- plyr::ddply(availability, .(year, gear), summarize,
                           lwr=quantile(availability, probs=.025),
                           upr=quantile(availability, probs=.975),
                           est=median(availability))
    availability2$space <- space
  } else {
    availability <- availability2 <- NULL
  }
  index.gear2$space <- index.strata$space <- space
  index.gear2$combinedoff <- index.strata$combinedoff <-
    availability2$combinedoff <- combinedoff
  index.gear2$fixlambda <- index.strata$fixlambda <-
    availability2$fixlambda <- fixlambda
  ## grab scenario from savedir
  scenario <- strsplit(savedir, split='/mcmc_')[[1]][2]
  out <- list(index.gear=index.gear2, index.strata=index.strata2,
              availability=availability2, scenario=scenario,
              R1_in=R1_in, R2_in=R2_in,
              PR1_in=PR1_in, PR2_in=PR2_in,
              PR1=PR1, PR2=PR2,
              beta1=beta1_tcn, beta2=beta2_tcn,
              R1_gcy=R1_gcy, R2_gcy=R2_gcy,
              D_gcyn=D_gcyn, covcor=covcor, savedir=savedir)
  saveRDS(out, file.path(savedir, 'results.mcmc.RDS'))
  return(out)
}


plot.posterior.predictive <- function(fit, results, plot=TRUE){
  ## Get some from each gear type and 0's and >0's
  message('Calculating posterior predictive...')
  x <- (1:nrow(Data_Geostat))
  dat <- Data_Geostat[, c("Gear", "Catch_KG")]
  R1 <- results$R1_in
  R2 <- results$R2_in
  ## Observation variances depend on the gear and sample
  if(model=='combined'){
    sigma.bts <- exp(as.data.frame(fit)[,'logSigmaM[1]']/1000)
    sigma.ats <- exp(as.data.frame(fit)[,'logSigmaM[2]']/1000)
  } else {
    sigma <- exp(as.data.frame(fit)[,'logSigmaM']/1000)
  }
  ## Genreate posterior predictive for each row of dat
  ppred <- array(NA, dim=c(nrow(R1), ncol(R1)))
  for(i in 1:nrow(dat)){
    ## Careful to use the right variance here
    if(model=='combined'){
      sig <- ifelse(dat$Gear[i]=='BT', sigma.bts, sigma.ats)
    } else {
      sig <- sigma
    }
    ppred[i,] <- exp(rnorm(n=ncol(R1), mean=log(R2[i,])-sig^2/2, sd=sig))*
      rbinom(n=ncol(R1), size=1, prob=R1[i,])
    ## convert to 0/1 for non encounters for easier plotting later
    if(dat$Catch_KG[i]==0)
      ppred[i,] <- ifelse(ppred[i,]==0, 0, 1)
  }
  message('Plotting posterior predictive...')
  ## Make plots of positive catches. Using percentile as metric
  Data_Geostat$percentile <- sapply(1:nrow(dat), function(i)
    mean(ppred[i,]<dat[i,'Catch_KG']))
  if(plot){
  for(zz in levels(Data_Geostat$Gear)){
    g <- ggplot(subset(Data_Geostat, Catch_KG>0 & Gear==zz),
                aes(x=Lon, y=Lat, color=percentile)) +
      scale_color_gradient2(midpoint=.5, low="blue", mid="white",
                            high="red", space ="Lab", limits=c(0,1) )+
      geom_jitter(height=.1, width=.1) + facet_wrap('Year') + theme_bw() +
      ## geom_tufteboxplot(median.type='line', hoffset=0, width=3)+
      ## facet_grid(Gear~year) +
      ## geom_point(aes(x=factor(rep), y=log(Catch_KG)), col='red', size=2) +
      ## xlab("Data number") + theme_bw() +
      ggtitle("Posterior predictive distribution for positive catches:", zz)
    ggsave(paste0(savedir, '/ppred_pos_', zz,'.png'), g, width=9, height=7)
  }
  ## Make some for the non-encounters. Better way to do this?
  ppred2 <- cbind(rep=1:nrow(dat), dat, year=Data_Geostat$Year, ppred) %>%
    gather(key=sample, value=catch, -rep, -Gear, -Catch_KG, -year)
  savedir <- results$savedir
  tmp <- ppred2 %>% filter(Catch_KG==0) %>% group_by(rep, Gear, year) %>%
    summarize(pct.zero=mean(catch==1))
  for(zz in levels(tmp$Gear)){
    g <- ggplot(subset(tmp, Gear==zz), aes(pct.zero))  + geom_histogram(bins=30) +
      facet_wrap('year') + xlab("Mean probability of encounter") +
      ggtitle("Posterior predictive distribution for non-encounters:", zz)
    ggsave(paste0(savedir, '/ppred_zeros_', zz,'.png'), g, width=9, height=5)
  }
  }
  return(invisible(Data_Geostat$percentile))
}

plot.covcor.mcmc <- function(results){
  ## Plot each one separately
  savedir <- results$savedir
  if(!is.null(results$covcor$covcor_omega1)){
    png(paste0(savedir, '/covcor_omega1.png'), width=7, height=5, res=500, units='in')
    plot.covcor(results$covcor$covcor_omega1, 'omega1')
    dev.off()
  }
  if(!is.null(results$covcor$covcor_omega2)){
    png(paste0(savedir, '/covcor_omega2.png'), width=7, height=5, res=500, units='in')
    plot.covcor(results$covcor$covcor_omega2, 'omega2')
    dev.off()
  }
  if(!is.null(results$covcor$covcor_epsilon1)){
    png(paste0(savedir, '/covcor_epsilon1.png'), width=7, height=5, res=500, units='in')
    plot.covcor(results$covcor$covcor_epsilon1, 'epsilon1')
    dev.off()
  }
}

plot.covcor <- function(covcor, Llab){
  par(mfrow=c(3,3), mgp=c(3,.5,0), mar=c(2,.5,2,.5), oma=c(0,0,2,0))
  for(i in 1:3){
    for(j in 1:3){
      if(j<=i){
        xlim <- c(-1,1)
        lab <- paste0('cor(',i, ',', j, ')')
        if(j==i){
          xlim <- c(0, max(sqrt(covcor[j,i,])))
          lab <- paste0('SD(', j,')')
          hist(sqrt(covcor[j,i,]), xlim=xlim, ylab=NA, xlab=NA, main=lab, yaxt='n'); box()
        } else {
          hist(covcor[j,i,], xlim=xlim, ylab=NA, xlab=NA, main=lab,
               yaxt='n'); box()
        }
      } else {
        plot(0,0, type='n', axes=FALSE, ann=FALSE)
      }
    }
  }
  mtext(Llab, line=0, outer=TRUE, cex=1.5)
}


plot.availability.map.mcmc <- function(results){
  savedir <- results$savedir
  if(is.null(results$D_gcyn)){
    warning("D_gcyn missing from index so skipping availability maps")
  } else {
    D_gcyn <- results$D_gcyn
    rm(results); gc(); gc() ## try to reduce memory usage
    Mapdetails <- make_map_info(Region, spatial_list=Spatial_List,
                                Extrapolation_List=Extrapolation_List)
    Mapdetails$Legend$x <- Mapdetails$Legend$x-70
    Mapdetails$Legend$y <- Mapdetails$Legend$y-45
    mdl <- Mapdetails
    Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
    Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
    ## For each draw calculate a surface of catchability by gear type
    MatTotal <- apply(D_gcyn, c(1,3,4), sum)
    ## Sum across first two strata
    MatBTS <- apply(D_gcyn[,-3,,], c(1,3,4), sum)/MatTotal
    ## Now the second two
    MatATS <- apply(D_gcyn[,-1,,], c(1,3,4), sum)/MatTotal
    ## Calculate CV and mean over posterior draws
    MatBTSCV <- apply(MatBTS, 1:2, function(x) sd(x)/mean(x))
    MatATSCV <- apply(MatATS, 1:2, function(x) sd(x)/mean(x))
    MatBTS <- apply(MatBTS, 1:2, mean)
    MatATS <- apply(MatATS, 1:2, mean)
    MatList <- list(BT=MatBTS, AT=MatATS)
    for(ii in c("AT", 'BT')){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatList[[ii]],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_availability_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=c(0,1),
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
    ## Repeat with CVs
    MatList <- list(BT=MatBTSCV, AT=MatATSCV)
    zlimtmp <- c(0, max(unlist(MatList)))
    for(ii in c("AT", 'BT')){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatList[[ii]],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_availability_CV_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=zlimtmp,
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
  }
}

plot.density.map.mcmc <- function(results){
  savedir <- results$savedir
  if(is.null(results$D_gcyn)){
    warning("D_gcyn missing from index so skipping density maps")
  } else {
    D_gcyn <- results$D_gcyn
    Mapdetails <- make_map_info(Region, spatial_list=Spatial_List,
                                Extrapolation_List=Extrapolation_List)
    Mapdetails$Legend$x <- Mapdetails$Legend$x-70
    Mapdetails$Legend$y <- Mapdetails$Legend$y-45
    mdl <- Mapdetails
    Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
    Years2Include = 1:length(Year_Set)#  which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
    ## For each strata calculate the mean log density
    MatStrata <- apply(log(D_gcyn), 1:3, mean)
    zlimtmp <- range(MatStrata)
    for(ii in 1:dim(MatStrata)[2]){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatStrata[,ii,],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_density_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=zlimtmp,
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
    ## Repeat with CVs
    MatStrata <- apply(D_gcyn, 1:3, function(x) sd(x)/mean(x))
    zlimtmp <- c(0, max(MatStrata))
    for(ii in 1:dim(MatStrata)[2]){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatStrata[,ii,],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_densityCV_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=zlimtmp,
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
  }
}


plot.R.map.mcmc <- function(results){
  Mapdetails <- make_map_info(Region, spatial_list=Spatial_List,
                              Extrapolation_List=Extrapolation_List)
  Mapdetails$Legend$x <- Mapdetails$Legend$x-70
  Mapdetails$Legend$y <- Mapdetails$Legend$y-45
  mdl <- Mapdetails
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = 1:length(Year_Set)
  if(is.null(results$R1_gcy)){
    warning("R1_gcy missing from index so skipping density maps")
  } else {
    ## For each strata calculate mean probability of occurence
    MatStrata <- results$R1_gcy
    zlimtmp <- c(0,1)
    for(ii in 1:dim(MatStrata)[[2]]){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatStrata[,ii,],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_R1_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=zlimtmp,
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
  }
  if(is.null(results$R2_gcy)){
    warning("R2_gcy missing from index so skipping density maps")
  } else {
    MatStrata <- log(results$R2_gcy)
    zlimtmp <- range(MatStrata)
    for(ii in 1:dim(MatStrata)[[2]]){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatStrata[,ii,],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/mcmc_map_R2_',ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=zlimtmp,
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Availability', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
  }
}

## Plot average pearson residuals by space
plot.pearson.mcmc <- function(results){
  dat <- cbind(Data_Geostat, PR1=results$PR1, PR2=results$PR2)
  for(zz in levels(dat$Gear)){
    g <- ggplot(subset(dat, Gear==zz & !is.na(PR2)),
                aes(Lon, Lat, size=abs(PR2), color=PR2<0)) +
      geom_point(alpha=.5) + facet_wrap('Year') + theme_bw() +
      ggtitle('Average Pearson resid for positive catches', zz)
    ggsave(paste0(savedir, '/pearson_pos_', zz, '.png'), g, width=9, height=7)
    ## Make sesne to look at binary ones?
    g <- ggplot(subset(dat, Gear==zz & !is.na(PR1)),
                aes(Lon, Lat, size=abs(PR1), color=PR1<0)) +
      geom_point(alpha=.5) + facet_wrap('Year') + theme_bw() +
      ggtitle('Average Pearson resid for non-encounters', zz)
    ggsave(paste0(savedir, '/pearson_enc_', zz, '.png'), g, width=9, height=7)
  }
  dat <- data.frame(year=Data_Geostat$Year, gear=Data_Geostat$Gear,
               PR2=results$PR2_in) %>% na.omit() %>%
    gather(rep, PR, -year, -gear) %>% mutate(pvalue=pnorm(PR))
  g <- ggplot(dat, aes(pvalue)) + geom_histogram(bins=50) +
    theme(axis.text.y=element_blank())
  ##  ggsave(paste0(savedir, '/pivotal_discrepancy.png'), g, width=9, height=4)
  ## g2 <- g+ facet_wrap('gear', scales='free')
  ## ggsave(paste0(savedir, '/pivotal_discrepancy_gear.png'), g2, width=9, height=4)
  g2 <- g+ facet_wrap(gear~year, scales='free')
  ggsave(paste0(savedir, '/pivotal_discrepancy_gear_year.png'), g2, width=9, height=7)
}


plot.betas.mcmc <- function(results, savedir){
  if(model !='combined'){
    dimnames(results$beta1) <- list(year=years, stratum=ifelse(model=='ats', strata.labels.ats, strata.labels.bts), iter=1:dim(results$beta1)[3])
    df1 <- cbind(beta='beta1', plyr::melt(results$beta1))
    dimnames(results$beta2) <- list(year=years, stratum=ifelse(model=='ats', strata.labels.ats, strata.labels.bts), iter=1:dim(results$beta2)[3])
    df2 <- cbind(beta='beta2', plyr::melt(results$beta2))
  } else {
    dimnames(results$beta1) <- list(year=years, stratum=strata.labels.combined, iter=1:dim(results$beta1)[3])
    df1 <- cbind(beta='beta1', plyr::melt(results$beta1))
    dimnames(results$beta2) <- list(year=years, stratum=strata.labels.combined, iter=1:dim(results$beta2)[3])
    df2 <- cbind(beta='beta2', plyr::melt(results$beta2))
  }
  df <- rbind(df1, df2) %>%
    plyr::ddply(.(stratum, beta, year), summarize,
          lwr=quantile(value, .025),
          upr=quantile(value, .975),
          med=median(value))
  g <-  ggplot(df, aes(year, med, fill=stratum, color=stratum)) +
  facet_grid(beta~stratum, scales='free_y')+
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
    geom_line(lwd=1.5) + geom_point()+ theme_bw() + ylab("value")
  ggsave(file.path(savedir, 'betas_mcmc.png'), g, width=7, height=5)
  ##   df$par.type <- sapply(strsplit(as.character(df$variable), split='\\['),
  ##                         function(x) x[[1]])
  ##   df$index <- as.numeric(sapply(strsplit(gsub('\\]', '', x=df$variable), split='\\['), function(x) x[[2]]))
  ##   temp <- data.frame(stratum=strata.labels.combined, year=rep(1:12, each=3), index=1:36)
  ##   df <- merge(df, temp, by='index')
  ##   df2 <- plyr::ddply(df, .(stratum, par.type, year), summarize,
  ##                lwr=quantile(value, .025),
  ##                upr=quantile(value, .975),
  ##                med=median(value))
  ##   g <-  ggplot(df2, aes(year, med, fill=stratum, color=stratum)) +  facet_grid(stratum~par.type)+
  ##     geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  ##     geom_line(lwd=1.5) + geom_point()+ theme_bw() + ylab("value")
  ##   ggsave(file.path(savedir, 'betas_mcmc.png'), g, width=7, height=5)
  ## }
}


plot.mcmc <- function(Obj, savedir, fit, n=8){
  results <- get.results.mcmc(Obj, fit)
  plot.index.mcmc(results, savedir)
  plot.betas.mcmc(results, savedir)
  plot.slow.mcmc(fit, savedir, n)
  plot.pairs.mcmc(fit, savedir)
  if(model=='combined')
    plot.availability.map.mcmc(results)
  plot.posterior.predictive(fit, results)
  ## plot.pearson.mcmc(results)
  plot.density.map.mcmc(results)
  plot.R.map.mcmc(results)
  plot.covcor.mcmc(results)
  ## Massage the output to get the beta's into a time format for ggplot
  pars.all <- names(fit)
  p <- pars.all[grep('lambda', x=pars.all)]
  if(length(p)>2){
    df <- plyr::melt(as.data.frame(fit)[,p], id.vars=NULL)
    df$par.type <- sapply(strsplit(as.character(df$variable), split='\\['),
                          function(x) x[[1]])
    df$index <- as.numeric(sapply(strsplit(gsub('\\]', '', x=df$variable), split='\\['), function(x) x[[2]]))
    temp <- data.frame(year=years, index=1:12)
    df <- merge(df, temp, by='index')
    df2 <- plyr::ddply(df, .(par.type, year), summarize,
                 lwr=quantile(value, .01),
                 upr=quantile(value, .99),
                 med=median(value))
    g <-  ggplot(df2, aes(year, med)) +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
      geom_line(lwd=2) + geom_point()+ theme_bw() + ylab("value")
    ggsave(file.path(savedir, 'lambdas_mcmc.png'), g, width=7, height=5)
  } else if(length(p)==1){
    ## just a constant but could be 1 or 2 parameters for the two LPs
    p <- pars.all[grep('lambda', x=pars.all)]
    df <- as.data.frame(fit)[,p, drop=FALSE] %>% gather(parameter, value)
    df2 <- df
    df2$parameter <- paste0('exp_', df2$parameter)
    df2$value <- exp(df2$value)
    df3 <- rbind(df, df2)
    g <- ggplot(df3, aes(value))  + geom_histogram(bins=20) +
      facet_wrap('parameter', scale='free') + theme_bw()
    ggsave(file.path(savedir, 'lambdas_mcmc.png'), g, width=7, height=5)
  } else {
    if(model=='combined') warning("not plot setup for 2 lambdas")
  }
  ## This is currently broken and probably not helpful anyway
  ## p <- pars.all[grep('Omegainput', x=pars.all)]
  ## if(length(p)>0){
  ##   df <- plyr::melt(as.data.frame(fit)[,p], id.vars=NULL)
  ##   df <- merge(df, temp, by='index')
  ##   n <- max(df$index)/2 # nmber of knots per factor (numer of rows)
  ##   temp <- data.frame(factor=rep(c('factor1', 'factor2'), each=n), knot=rep(1:n,times=2), index=1:(2*n))
  ##   df2 <- plyr::ddply(df, .(factor, par.type, knot), summarize,
  ##                lwr=quantile(value, .01),
  ##                upr=quantile(value, .99),
  ##                med=median(value))
  ##   g <-  ggplot(df2, aes(knot, med, fill=factor, color=factor)) +  facet_grid(factor~par.type)+
  ##     geom_pointrange(aes(ymin=lwr, ymax=upr), alpha=.5) +
  ##     geom_point()+ theme_bw() + ylab("value")
  ##   ggsave(file.path(savedir, 'omegas_mcmc.png'), g, width=7, height=5)
  ## }
}





plot.index.mcmc <- function(results, savedir){
  if(is.null(results)){ message("index is NULL so skipping plots"); return()}
  g <- ggplot(results$index.gear, aes(year, y=est, color=gear, group=gear, fill=gear)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
    geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
    ylab('log abundance')
  ggsave(file.path(savedir, 'index_gear_mcmc.png'), g, width=7, height=5)
  g <- ggplot(results$index.strata, aes(year, y=est, color=stratum,  group=stratum, fill=stratum)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
    geom_line(lwd=1.5, alpha=.5) + theme_bw() + # facet_wrap('stratum')+
    ylab('log abundance')
  ggsave(file.path(savedir, 'index_strata_mcmc.png'), g, width=7, height=5)
  if(model=='combined'){
    g <- ggplot(results$availability, aes(year, y=est, color=gear, group=gear, fill=gear)) +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
      geom_line(lwd=1.5, alpha=.5) + theme_bw() +# facet_wrap('gear')+
      ylab('Availability to gear') + ylim(0,1)
    ggsave(file.path(savedir, 'availability_mcmc.png'), g, width=7, height=5)
    ## Do relative densities by strata and year for median
    tmp <- dcast(results$index.strata, year~stratum, value.var='est')
    tmp[,2:4] <- exp(tmp[,2:4])/rowSums(exp(tmp[2:4]))
    index.strata.pct <- plyr::melt(tmp, 'year', variable.name='stratum',
                             value.name='pct.density')
    index.strata.pct$stratum <- factor(index.strata.pct$stratum,
                                       levels=rev(levels(results$index.strata$stratum)) )
    g <- ggplot(index.strata.pct, aes(year, pct.density, fill=stratum)) +
      geom_area()
    ggsave(file.path(savedir, 'pct_strata_mcmc.png'), g, width=7, height=5)
  }
}

plot.slow.mcmc <- function(fit, savedir, n=8){
  mon <- monitor(fit, print=FALSE)
  mon <- as.data.frame(mon)
  mon$par <- row.names(mon)
  mon.summary <- summarize(mon, minESS=min(n_eff), maxRhat=max(Rhat),
                           minESSBulk=min(Bulk_ESS),
                           minESSTail=min(Tail_ESS))
  mon.summary$ndivs <- get_num_divergent(fit)
  mon.summary$nmaxtd <- get_num_max_treedepth(fit)
  write.csv(x=mon.summary, file='monitor_summary.csv')
  mon$par.type <- 'fixed'
  mon$par.type[grep('Omegainput|Epsiloninput', mon$par)] <- 'random'
  mon$par.name <- sapply(strsplit(mon$par, split='\\['), function(x) x[[1]])
  mon$space <- space; mon$combinedoff <- combinedoff
  mon$fixlambda <- fixlambda
  mon <- mon[order(mon$n_eff),]
  ## row.names(mon) <- NULL;
  mon$energy__ <- NULL
  pars.slow.fixed <- mon[mon$par.type=='fixed','par'][1:n]
  pars.slow.random <- mon[mon$par.type=='random','par'][1:n]
  print(mon[1:10,c('n_eff', 'Rhat')])
  g <- ggplot(mon, aes(x=0, y=n_eff, color=par.name)) + geom_jitter(alpha=.5) +
    scale_y_log10() + facet_wrap('par.type') + ylab('log10 ESS') +
    theme(axis.text.x=element_blank())
  ggsave(file.path(savedir, 'ess.png'), g, width=7, height=5)
  saveRDS(file.path(savedir, 'monitor.RDS'), object=mon)
  png(paste0(savedir, '/pairs_slow_fixed.png'), width=7, height=5, res=500, units='in')
  pairs(fit, pars=pars.slow.fixed, gap=0)
  dev.off()
  if(length(na.omit(pars.slow.random))>2){
    png(paste0(savedir, '/pairs_slow_random.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=pars.slow.random, gap=0)
    dev.off()
  }
}

plot.pairs.mcmc <- function(fit, savedir){
  message("Making pairs plots which can be slow..")
  pars.all <- names(fit)
  p <- pars.all[grep('L_omega1_z', x=pars.all)]
  if(length(p)>1){
    png(paste0(savedir, '/pairs_L_omega1.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=p, gap=0)
    dev.off()
  }
  p <- pars.all[grep('L_omega2_z', x=pars.all)]
  if(length(p)>1){
    png(paste0(savedir, '/pairs_L_omega2.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=p, gap=0)
    dev.off()
  }
  p <- pars.all[grep('L_epsilon1_z', x=pars.all)]
  if(length(p)>1){
    png(paste0(savedir, '/pairs_L_epsilon1.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=p, gap=0)
    dev.off()
  }
  p <- pars.all[grep('L_epsilon2_z', x=pars.all)]
  if(length(p)>1){
    png(paste0(savedir, '/pairs_L_epsilon2.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=p, gap=0)
    dev.off()
  }
  p <- pars.all[grep('kappa|Sigma|lp__|lambda|rho', x=pars.all)]
  png(paste0(savedir, '/pairs_fixed.png'), width=7, height=5, res=500, units='in')
  pairs(fit, pars=p, gap=0)
  dev.off()
  p <- pars.all[grep('lambda|Beta_mean|L_beta', x=pars.all)]
  if(length(p)<15 & length(p) >0){
    png(paste0(savedir, '/pairs_scale.png'), width=7, height=5, res=500, units='in')
    pairs(fit, pars=p, gap=0)
    dev.off()
  } else {
    warning("in plot.pairs.mcmc the 'scale' plot had too many/few parameters")
  }
}

calculate.index <- function(Opt, Report, model, space, log, strata){
  ## If available use the bias corrected versions
  if(is.null(Opt$SD)){
    ## If no sdeport need to exit with just the MLEs
    if(log){
      if(strata){
        ests <- t(log(Report$Index_cyl[,,1]))
      } else {
        ests <- t(Report$ln_ColeIndex_cy)
      }
    } else {
      if(strata){
        ests <- t(Report$Index_cyl)
      } else {
        ests <- t(Report$ColeIndex_cy)
      }
    }
    ests <- as.numeric(ests)
    ses <- rep(NA, len=length(ests))
  } else {
    ## SD exists
    if(!is.null(Opt$SD$unbiased)){
      rr <- as.list(Opt$SD, "Est. (bias.correct)", report=TRUE)
    } else {
      ## Not available so used uncorrected version
      rr <- as.list(Opt$SD, what='Estimate', report=TRUE)
    }
    ## Assuming for now that the SDs aren't bias corrected
    ss <- as.list(Opt$SD, what='Std. Error', report=TRUE)
    if(strata){
      ## the strata versions
      if(log){
        ## calculate in log space?
        ses <- t(ss$ln_Index_cyl[,,1])
        ests <- t(rr$ln_Index_cyl[,,1])
      } else {
        ses <- t(ss$Index_cyl[,,1])
        ests <- t(rr$Index_cyl[,,1])
      }
      snames <- c('stratum1', 'stratum2', 'stratum3')
    } else {
      ## the versions the gear sees
      if(log){
        ## calculate in log space?
        ses <- t(ss$ln_ColeIndex_cy)
        ests <- t(rr$ln_ColeIndex_cy)
      } else {
        ses <- t(ss$ColeIndex_cy)
        ests <- t(rr$ColeIndex_cy)
      }
      snames <- c('total', 'bts', 'ats')
    }
  }
  ## Chop of years of missing ATS if necessary
  yrs <- years[which(min(years):max(years) %in% years)]
  if(model=='combined'){
    snames <- c('total', 'bts', 'ats')
    if(strata)
      snames <- c('stratum1', 'stratum2', 'stratum3')
    index <- data.frame(model=model, space=space,  year=yrs,
                        strata=rep(snames, each=length(years)),
                        est=as.vector(ests), se=as.vector(ses))
  } else {
    ## ATS or BTS is just a single column
    index <- data.frame(model=model, space=space,  year=yrs,
                        strata=model,
                        est=as.vector(ests), se=as.vector(ses))
  }
  index <- within(index, {lwr <- est-1.96*se; upr <- est+1.96*se})
  return(index)
}

calculate.index.old <- function(Opt, Report, model, space, log, strata){
  if(log){
    ## calculate in log space?
    if(!strata)
      stop("Doesnt make sense to have combined model in log space")
    tmp <- which(names(Opt$SD$value) %in% 'ln_Index_cyl')
    ii <- log(Report$Index_cyl[,,1])
  } else {
    tmp <- which(names(Opt$SD$value) %in% 'Index_cyl')
    ii <- Report$Index_cyl[,,1]
  }
  index <- data.frame(model=model, space=space,  year=years)
  if(model=='combined'){
    if(!strata){
      ## Manually calculate SE for the total biomass index. Since it's a sum of
      ## the three the derivatives are all 1 and so the SE is the sqrt of the sum
      ## of all of the variances and covariances. This feature is not coded
      ## into VAST yet so have to do it manually. also note that the order of
      ## the Index_cyl matrix in vector form is Index_11, Index_21, Index_31,
      ## Index_12,.. etc. This effects the subsetting below
      cov.index <- Opt$SD$cov[tmp,tmp]
      ## combined is all three strata summed
      index1 <- data.frame(index, strata='total',  est=apply(ii[1:3,], 2, sum),
                           se=sqrt(sapply(1:nyr, function(i) {j=1:3+3*(i-1);
                             sum(cov.index[j,j])})))
      ## sum the top two to get what the ATS sees
      index2 <- data.frame(index, strata='ats', est=apply(ii[2:3,], 2, sum),
                           se=sqrt(sapply(1:nyr, function(i) {j=(1:3+3*(i-1))[-1];
                             sum(cov.index[j,j])})))
      ## likewise the BTS is just the first strata
      index3 <- data.frame(index, strata='bts', est=apply(ii[1:2,], 2, sum),
                           se=sqrt(sapply(1:nyr, function(i) {j=(1:3+3*(i-1))[-3];
                             sum(cov.index[j,j])})))
      index <- rbind(index1,index2, index3)
    } else {
      ## If we just want to get the 3 strata without the summation
      sdtmp <- matrix(sqrt(diag(Opt$SD$cov[tmp,tmp])), nrow=3)
      index1 <- data.frame(index, strata='strata1', est=ii[1,], se=sdtmp[1,])
      index2 <- data.frame(index, strata='strata2', est=ii[2,], se=sdtmp[2,])
      index3 <- data.frame(index, strata='strata3', est=ii[3,], se=sdtmp[3,])
      index <- rbind(index1,index2, index3)
    }
  } else {
    ## Either BTS or ATS so single
    ## chop off missing years for ATS case
    tmp2 <- which(min(years):max(years) %in% years)
    tmp <- tmp[tmp2]
    index <- data.frame(index, strata=model, est=Opt$SD$value[tmp], se=Opt$SD$sd[tmp])
  }
  index <- within(index, {lwr <- est-1.96*se; upr <- est+1.96*se})
  return(index)
}


plot.vastfit <- function(results, savedir,  plotQQ=FALSE, plotmaps=FALSE){
  beta1 <- results$Report$beta1_tc
  beta2 <- results$Report$beta2_tc
  if(model !='combined'){
    dimnames(beta1) <- list(year=years, stratum=ifelse(model=='ats', strata.labels.ats, strata.labels.bts))
    df1 <- cbind(beta='beta1', plyr::melt(beta1))
    dimnames(beta2) <- list(year=years, stratum=ifelse(model=='ats', strata.labels.ats, strata.labels.bts))
    df2 <- cbind(beta='beta2', plyr::melt(beta2))
  } else {
    dimnames(beta1) <- list(year=years, stratum=strata.labels.combined)
    df1 <- cbind(beta='beta1', plyr::melt(beta1))
    dimnames(beta2) <- list(year=years, stratum=strata.labels.combined)
    df2 <- cbind(beta='beta2', plyr::melt(beta2))
  }
  betas <- rbind(df1,df2)
  if(nrow(betas)>0){
    g <-  ggplot(betas, aes(year, value, group=beta, color=beta)) +
      facet_wrap(stratum~.)+ geom_line(lwd=2)
    ggsave(file.path(savedir, 'betas.png'), g, width=7, height=5)
  }

  ## Need to reconstruct the Density including the log-normal bias
  ## adjustment
  lambdas <- subset(results$est, par=='lambda1_k')
  if(nrow(lambdas)>3){
    lambdas$year <- years
    g <- ggplot(lambdas, aes(year, y=est)) +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
      geom_line() + geom_point()+ theme_bw() + ylab('BTS catchability')
    ggsave(file.path(savedir, 'tv_lambdas.png'), g, width=7, height=5)
  }
  sigtmp <- results$Report$SigmaM[as.numeric(Data_Geostat$Gear)]^2/2
  df <- data.frame(obs=log(Data_Geostat$Catch_KG),
                   predicted= log(results$Report$R2_i)-sigtmp,
                   gear=Data_Geostat$Gear,
                   year=Data_Geostat$Year)
  df <- subset(df, obs>0) ## drop zeroes
  g <- ggplot(df, aes(obs, predicted)) + facet_grid(gear~year) + geom_point(alpha=.5) +
    geom_abline(slope=1, intercept=0)
  ggsave(file.path(savedir, 'obs_vs_pred.png'), g, width=10, height=5)
  ## ## This was old code to look at the input fields but are not really
  ## ## helpful anymore
  ## if(results$Index$space[1]!="NS"){
  ##   fields <- data.frame(model=results$Index$model[1], space=results$Index$space[1],
  ##                        omegainput1=results$Report$Omegainput1_sf,
  ##                        omega1=results$Report$Omega1_sc,
  ##                        omegainput2=results$Report$Omegainput2_sf,
  ##                        omega2=results$Report$Omega2_sc,
  ##                        E_km=results$Inputs$loc$E_km,
  ##                        N_km=results$Inputs$loc$N_km)
  ##   fields.long <- plyr::melt(fields, id.vars=c('model', 'space', 'E_km', 'N_km'),
  ##                       factorsAsStrings=FALSE)
  ##   if(results$Index$model[1]=='combined'){
  ##     fields.long$strata <- paste0('strata_',unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
  ##                                                          function(x) x[2])))
  ##   } else {
  ##     fields.long$strata <- results$Index$model[1]
  ##   }
  ##   fields.long$type <- unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
  ##                                     function(x) x[1]))
  ##   fields.long$component <- 'Component=1'
  ##   fields.long$component[grep(fields.long$type, pattern='2')] <- 'Component=2'
  ##   fields.long$type <- gsub("1|2", "", x=fields.long$type)
  ##   fields.long$type <- factor(fields.long$type, levels=c('omegainput', 'omega'))
  ##   fields.long <- plyr::ddply(fields.long, .(type, space, component), mutate,
  ##                        normalized=value/sd(value))
  ##   Col  <-  colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
  ##   g <- ggplot(fields.long, aes(E_km, N_km, col=value)) +
  ##     geom_point(size=1) +
  ##     facet_grid(component+type~strata) +
  ##     scale_colour_gradientn(colours = Col(15)) + theme_bw()
  ##   ggsave(file.path(savedir, 'map_omegas.png'), g, width=9, height=6, units='in')
  ## }
  Report <- results$Report
  g <- ggplot(results$Index, aes(year, y=est, group=strata,
                                 fill=strata), color=strata) +
    geom_ribbon(aes(ymin=est-1.96*se, ymax=est+1.96*se), alpha=.5) +
    geom_line() + geom_point()+ theme_bw() +
    ylab('log abundance')
  ggsave(file.path(savedir, 'index.png'), g, width=7, height=5)
  ## Also create an index of the individual strata
  if(results$model=='combined'){
    g <- ggplot(results$Index.strata,
                aes(year, y=est, group=strata, fill=strata, color=strata)) +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
      geom_line() + geom_point()+ theme_bw() + ylab('log abundance')
    ggsave(file.path(savedir, 'index_strata.png'), g, width=7, height=5)
  ## Do relative densities by strata and year for median
  tmp <- dcast(results$Index.strata, year~strata, value.var='est')
  tmp[,2:4] <- exp(tmp[,2:4])/rowSums(exp(tmp[2:4]))
  index.strata.pct <- plyr::melt(tmp, 'year', variable.name='strata',
                           value.name='pct.density')
  index.strata.pct$strata <- factor(index.strata.pct$strata,
                                     levels=rev(levels(results$Index.strata$strata)) )
  g <- ggplot(index.strata.pct, aes(year, pct.density, fill=strata)) +
    geom_area()
  ggsave(file.path(savedir, 'pct_strata.png'), g, width=7, height=5)
  }
  Mapdetails <- make_map_info(Region, spatial_list=Spatial_List,
                              Extrapolation_List=Extrapolation_List)
  Mapdetails$Legend$x <- Mapdetails$Legend$x-70
  Mapdetails$Legend$y <- Mapdetails$Legend$y-45
  mdl <- Mapdetails
  ## This was causing problems and not sure why. Will fix later.
  if(TmbData$n_c>1){
    Report$D_xcy <- Report$D_gcy
    plot_factors(Report, ParHat=results$ParHatList, Data=TmbData, SD=Opt$SD,
                 mapdetails_list=Mapdetails, plotdir=savedir)
  }
  Enc_prob <- plot_encounter_diagnostic(Report=Report,
                                        Data_Geostat=Data_Geostat,
                                        DirName=savedir)
  ## Decide which years to plot
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  TmbData$n_x <- TmbData$n_g
  if(plotQQ){
    Q <- plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                                  FileName_Phist="Posterior_Predictive-Histogram",
                                  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=savedir )
    ## MapDetails_List = make_map_info( "Region"=Region,
    ##                                 "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
    ##                                 "Extrapolation_List"=Extrapolation_List )
    plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'],
                   TmbData=TmbData, Report=Report, Q=Q, savedir=savedir,
                   MappingDetails=Mapdetails[["MappingDetails"]],
                   PlotDF=Mapdetails[["PlotDF"]],
                   MapSizeRatio=Mapdetails[["MapSizeRatio"]],
                   Xlim=Mapdetails[["Xlim"]],
                   Ylim=Mapdetails[["Ylim"]], FileName=savedir,
                   Year_Set=Year_Set, Years2Include=Years2Include,
                   Rotate=Mapdetails[["Rotate"]],
                   Cex=Mapdetails[["Cex"]],
                   Legend=Mapdetails[["Legend"]],
                   zone=Mapdetails[["Zone"]], mar=c(0,0,2,0),
                   oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  }
  if(plotmaps){
    plot_anisotropy( FileName=paste0(savedir,"/Aniso.png"), Report=Report,
                    TmbData=TmbData )
    ## Some built-in maps
    tmp <- c(1,2,3, 11, 12)
    if(results$Index$space[1]=='ST') tmp <- c(tmp, 6,7)
    ## Temporary hack to get v8.0.0 to work with this function
    Report$D_xcy <- Report$D_gcy
    Report$R1_xcy <- Report$R1_gcy
    Report$R2_xcy <- Report$R2_gcy
    Report$D_gcy <- Report$R1_gcy <- Report$R2_gcy <- NULL
    TmbData$X_xtp <- TmbData$X_gtp
    Dens_xt = plot_maps(plot_set=tmp,
                        ## MappingDetails=Mapdetails[["MappingDetails"]],
                        Report=Report, Sdreport=Opt$SD,
                        TmbData=TmbData,
                        PlotDF=Mapdetails[["PlotDF"]],
                        MapSizeRatio=Mapdetails[["MapSizeRatio"]],
                        ## Xlim=Mapdetails[["Xlim"]],
                        ## Ylim=Mapdetails[["Ylim"]],
                        working_dir=paste0(savedir,'/'),
                        Year_Set=Year_Set, Years2Include=Years2Include,
                        ## Rotate=Mapdetails[["Rotate"]],
                        ## Cex=Mapdetails[["Cex"]],
                        ##Legend=Mapdetails[["Legend"]],
                        ## zone=Mapdetails[["Zone"]],
                        mar=c(0,0,2,0),
                        oma=c(3.5,3.5,0,0), cex=1.8)
    ## Dens_DF = cbind( "Density"=as.vector(Dens_xt),
    ##                 "Year"=Year_Set[col(Dens_xt)],
    ##                 "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'],
    ##                 "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )
    if(!is.null(results$Opt$SD)){
      Index = plot_biomass_index( DirName=savedir, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
      ##  pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )
      plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=savedir, Year_Set=Year_Set)
    }
    if(results$Index$model[1]=='combined'){
      ## Plot ratio of observed/predicted by grid cell for the three gear types
      MatDat <- (tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                        FUN=mean, na.rm=TRUE))
      MatExp <- Report$D_xcy
      MatRatio <- array(NA, dim=dim(MatDat))
      ## BTS is sum over first two strata
      MatRatio[,1,] <- MatDat[,1,]/apply(MatExp[,-3,], 2, sum)
      MatRatio[,2,] <- MatDat[,2,]/MatExp[,2,]
      MatRatio[,3,] <- MatDat[,3,]/MatExp[,3,]
      MatRatio <- log(MatRatio)
      MatRatio[is.infinite(MatRatio)] <- NA
      for(ii in 1:3){
        PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                   Mat=MatRatio[,ii,Years2Include,drop=TRUE],
                   PlotDF=mdl$PlotDF,
                   MapSizeRatio=mdl$MapSizeRatio,
                   Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                   FileName=paste0(savedir, 'map_data_ratio_', ii),
                   Year_Set=Year_Set[Years2Include],
                   Legend=mdl$Legend, zlim=range(MatRatio, na.rm=TRUE),
                   mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                   textmargin='log(Obs/Exp)', zone=mdl$Zone, mar=c(0,0,2,0),
                   oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
      }
    }
  }
  ## Pearson resids for detection and catch rate
  D_i <- Report$R1_i*Report$R2_i
  PR1_i <- PR2_i <- rep(NA, length(D_i))
  for(i in 1:length(D_i)){
    ## bernoulli for presence
    mui <- Report$R1_i[i]
    obs <- as.numeric(Data_Geostat$Catch_KG[i]>0)
    PR1_i[i] <- (obs-mui)/sqrt(mui*(1-mui)/1)
    ## log-normal for catch rate; NA for 0 observations
    obs <- Data_Geostat$Catch_KG[i]
    if(obs>0){
      ## make sure to use the right variance as this depends on gear type
      gr <- as.numeric(Data_Geostat$Gear[i])
      PR2_i[i] <- (log(obs)-log(Report$R2_i[i])+Report$SigmaM[gr]^2/2)/Report$SigmaM[gr]
    }
  }
  df <- cbind(Data_Geostat, PR1=PR1_i, PR2=PR2_i, positive=ifelse(Data_Geostat$Catch_KG>0,1,0))
  xlim <- range(df$Lon); ylim <- range(df$Lat)
  tmp <- 1:3
  if(results$model!='combined') tmp <- 1
  for(gr in tmp){
    gt <- levels(Data_Geostat$Gear)[gr]
    g <- ggplot(subset(df, Gear==gt & positive==1), aes(Lon, Lat, size=abs(PR2), color=PR2>0))+
      geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
      scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
    ggsave(filename=paste0(savedir, '/Pearson_resid_catchrate_', gr, '.png'), plot=g,
           width=7, height=5)
    g <- ggplot(subset(df, Gear==gt), aes(Lon, Lat, size=abs(PR1), color=PR1>0))+
      geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
      scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
    ggsave(filename=paste0(savedir, '/Pearson_resid_encounter_', gr, '.png'), plot=g,
           width=7, height=5)
  }


  ## QQplots
  get.qq <- function(gr, presence){
    gears <- levels(Data_Geostat$Gear)
    qq <- ldply(years, function(i){
      tmp <- subset(df, Year==i & positive==1 & Gear==gears[gr])
      if(nrow(tmp)>0) {
        qq <- data.frame(qqnorm(tmp$PR2, plot.it=FALSE))
        x <- data.frame(year=i, gear=gears[gr], qq)
        return(x)
      }
    })
  }
  if(results$model=='combined'){
    qq <- rbind(get.qq(1,1), get.qq(2,1), get.qq(3,1))
  } else {
    qq <- get.qq(1,1)
  }
  g <- ggplot(qq, aes(x,y, group=gear, color=gear)) +
    geom_abline(slope=1,intercept=0) + facet_wrap('year') +
    geom_point(alpha=.5) + theme_bw() + xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles")
  ggsave(filename=paste0(savedir, '/QQplot_catchrate.png'), plot=g,
         width=7, height=5)
  ## if(results$model=='combined')
    ## plot.change(Report)
}

##   ## This is a modified version of plot_residuals meant to work with my
## ## combined strata model

## observed.catches.by.geartype <- function(gear){
##   obs_rate_xy <- total_num_xy <- obs_num_xy <-
##     matrix(NA, nrow=TmbData$n_x, ncol=nrow(TmbData$t_yz) )
##   ## Loop through each year and calculate obs for each cell. Be careful to
##   ## only select the appropriate gear type and then sum across the correct strata
##   for(yr in  1:nrow(TmbData$t_yz)){
##     tmp <- subset(Data_Geostat, Year== years[yr] & Gear == gear)
##     if(nrow(tmp)>0){
##       index.tmp <- factor(tmp$knot_i, levels=1:TmbData$n_x-1)
##       obs_rate_xy[,yr] <- tapply(tmp$Catch_KG>0, INDEX=index.tmp, FUN=mean )
##       total_num_xy[,yr] <- tapply(tmp$Catch_KG, INDEX=index.tmp, FUN=length)
##     } else {
##       total_num_xy[,yr] <- NA
##     }
##   }
##   return(list(obs_rate_xy=obs_rate_xy, total_num_xy=total_num_xy))
## }

## Q1_xy <- list()
## #### First get the Pearon resids for BT
## ## Sum across P1 for the first two strata then calculate R1 manually. I
## ## don't think there's another way to do this.
## exp_rate_xy <- 1-exp(-exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)))
## temp <- observed.catches.by.geartype('BT')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q1_xy[[1]] <- (obs_num_xy-exp_num_xy)/sqrt(exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS1
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('AT2')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q1_xy[[2]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS2
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('AT3')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q1_xy[[3]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )

## textmargin = "Pearson residual"
## Col = colorRampPalette(colors=c("blue","white","red"))
## for( zI in 1:3 ){
##   ## Sometimes expected is 1 to machine precision which causes infinite
##   ## residual. For now truncating to 3
##   ## for(ii in 1:3) Q1_xy[[ii]][which(is.infinite(Q1_xy[[ii]]), arr.ind=TRUE)] <- NA
##   ## for(ii in 1:3) Q1_xy[[ii]][which(abs(Q1_xy[[ii]])>50, arr.ind=TRUE)] <- NA
##   Q1_xy[[zI]]  <-  ifelse( abs(Q1_xy[[zI]])>3, 3*sign(Q1_xy[[zI]]), Q1_xy[[zI]] )
##   zlim  <-  c(-1,1) *3#ceiling(max(abs(unlist(Q1_xy)),na.rm=TRUE))
##   PlotMap_Fn( MappingDetails=mdl$MappingDetails, Mat=Q1_xy[[zI]],
##              PlotDF=mdl$PlotDF, Col=Col, zlim=zlim, ignore.na=TRUE,
##              MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
##              FileName=paste0(savedir, '/Pearson_resid_presence_', zI),
##              Year_Set=paste0("Residuals--",1:ncol(Q1_xy[[zI]])),
##              zone=mdl$Zone, Legend=mdl$Legend,
##              mfrow=c(ceiling(sqrt(length(Years2Include))),
##                      ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
##              mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
## }


## Q2_xy <- list()
## #### First get the Pearon resids for BT
## ## Sum across P1 for the first two strata then calculate R1 manually. I
## ## don't think there's another way to do this.
## r1temp <- 1-exp(-exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)))
## ## For R2=exp(p1+p2)/r1
## exp_rate_xy <-
##   exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)+apply(Report$P2_xcy[,-3,], c(1,3), sum))/r1temp
## temp <- observed.catches.by.geartype('BT')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q2_xy[[1]] <- (obs_num_xy-exp_num_xy)/sqrt(exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS1
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('AT2')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q2_xy[[2]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS2
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('AT3')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q2_xy[[3]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )

## for( zI in 1:3 ){
##   ## Sometimes expected is 1 to machine precision which causes infinite
##   ## residual. For now truncating to 3
##   ## for(ii in 1:3) Q2_xy[[ii]][which(is.infinite(Q2_xy[[ii]]), arr.ind=TRUE)] <- NA
##   ## for(ii in 1:3) Q2_xy[[ii]][which(abs(Q2_xy[[ii]])>50, arr.ind=TRUE)] <- NA
##   Q2_xy[[zI]]  <-  ifelse( abs(Q2_xy[[zI]])>3, 3*sign(Q2_xy[[zI]]), Q2_xy[[zI]] )
##   zlim  <-  c(-1,1) *3#ceiling(max(abs(unlist(Q2_xy)),na.rm=TRUE))
##   PlotMap_Fn( MappingDetails=mdl$MappingDetails, Mat=Q2_xy[[zI]],
##              PlotDF=mdl$PlotDF, Col=Col, zlim=zlim, ignore.na=TRUE,
##              MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
##              FileName=paste0(savedir, '/Pearson_resid_presence_', zI),
##              Year_Set=paste0("Residuals--",1:ncol(Q2_xy[[zI]])),
##              zone=mdl$Zone, Legend=mdl$Legend,
##              mfrow=c(ceiling(sqrt(length(Years2Include))),
##                      ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
##              mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
## }



## ### Method #3 -- Pearson residuals
##   sum_obs_xy = sum_exp_xy = var_exp_xy = matrix(NA, nrow=TmbData$n_x, ncol=nrow(TmbData$t_yz) )
##   for( yI in 1:nrow(TmbData$t_yz) ){
##     which_i_in_y = ( TmbData$t_iz == outer(rep(1,TmbData$n_i),TmbData$t_yz[yI,]) )
##     which_i_in_y = which( apply(which_i_in_y,MARGIN=1,FUN=all) )
##     which_i_in_y_and_pos = intersect( which_i_in_y, which_pos )
##     which_ipos_in_y = ( TmbData$t_iz[which_pos,] == outer(rep(1,length(which_pos)),TmbData$t_yz[yI,]) )
##     which_ipos_in_y = which( apply(which_ipos_in_y,MARGIN=1,FUN=all) )
##     if( length(which_i_in_y_and_pos)>0 ){
##       sum_obs_xy[,yI] = tapply( TmbData$b_i[which_i_in_y_and_pos], INDEX=factor(TmbData$s_i[which_i_in_y_and_pos],levels=1:TmbData$n_x-1), FUN=sum )
##       sum_exp_xy[,yI] = tapply( bpred_ipos[which_ipos_in_y], INDEX=factor(TmbData$s_i[which_i_in_y_and_pos],levels=1:TmbData$n_x-1), FUN=sum )
##       var_exp_xy[,yI] = tapply( bvar_ipos[which_ipos_in_y], INDEX=factor(TmbData$s_i[which_i_in_y_and_pos],levels=1:TmbData$n_x-1), FUN=sum )
##     }
##   }
##   Q2_xy = (sum_obs_xy - sum_exp_xy) / sqrt(var_exp_xy)





