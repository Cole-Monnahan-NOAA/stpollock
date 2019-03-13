

library(TMB)
## library(devtools)
## remove.packages('VAST')
## install_github('James-Thorson/VAST')
library(VAST)
library(reshape2)
library(ggplot2)
library(plyr)
library(TMBhelper)
library(snowfall)
library(maps)
library(mapdata)
library(abind)
Version <- "VAST_v7_0_0"

source("simulator.R")


process.results <- function(Opt, Obj, Inputs, model, space, savedir){
  Report  <-  Obj$report()
  ParHatList <- Obj$env$parList(Opt$par)
  ParHat <- Opt$par
  if(is.null(Opt$SD$cov.fixed)){ # no SD
    SE <- rep(Inf, len=length(ParHat))
  } else {
    SE <- sqrt(diag(Opt$SD$cov.fixed))
  }
  est <- data.frame(par=names(ParHat), est=ParHat, lwr=ParHat-1.96*SE,
                    upr=ParHat+1.96*SE)
  est$significant <- !(est$lwr<0 & est$upr>0)
  write.csv(est, file=paste0(savedir, "/estimates.csv"))
  Index <- calculate.index(Opt, Report, model, space, log=FALSE, strata=FALSE)
  Index.strata <- calculate.index(Opt, Report, model, space, log=TRUE, strata=TRUE)
  Save  <-  list(Index=Index, Opt=Opt, Report=Report, ParHat=ParHat,
                 ParHatList=ParHatList, est=est, Index.strata=Index.strata,
                 SE=SE, Inputs=Inputs, savedir=savedir)
  save(Save, file=paste0(savedir,"/Save.RData"))
  return(Save)
}

calculate.index <- function(Opt, Report, model, space, log, strata){
  ## If available use the bias corrected versions
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
  ## Chop of years of missing ATS if necessary
  yrs <- years[which(min(years):max(years) %in% years)]
  if(model=='combined'){
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

plot.vastfit <- function(results){
  df <- data.frame(obs=Data_Geostat$Catch_KG,
                   predicted=results$Report$R2_i, gear=Data_Geostat$Gear,
                   year=Data_Geostat$Year)
  df <- subset(df, obs>0) ## drop zeroes
  g <- ggplot(df, aes(log(obs), log(predicted))) + facet_grid(gear~year) + geom_point(alpha=.5) +
    geom_abline(slope=1, intercept=0)
  ggsave(file.path(savedir, 'obs_vs_pred.png'), g, width=10, height=5)
  if(results$Index$space[1]!="NS"){
    fields <- data.frame(model=results$Index$model[1], space=results$Index$space[1],
                         omegainput1=results$Report$Omegainput1_sf,
                         omega1=results$Report$Omega1_sc,
                         omegainput2=results$Report$Omegainput2_sf,
                         omega2=results$Report$Omega2_sc,
                         E_km=results$Inputs$loc$E_km,
                         N_km=results$Inputs$loc$N_km)
    fields.long <- melt(fields, id.vars=c('model', 'space', 'E_km', 'N_km'),
                        factorsAsStrings=FALSE)
    if(results$Index$model[1]=='combined'){
      fields.long$strata <- paste0('strata_',unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
                                                           function(x) x[2])))
    } else {
      fields.long$strata <- results$Index$model[1]
    }
    fields.long$type <- unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
                                      function(x) x[1]))
    fields.long$component <- 'Component=1'
    fields.long$component[grep(fields.long$type, pattern='2')] <- 'Component=2'
    fields.long$type <- gsub("1|2", "", x=fields.long$type)
    fields.long$type <- factor(fields.long$type, levels=c('omegainput', 'omega'))
    fields.long <- ddply(fields.long, .(type, space, component), mutate,
                         normalized=value/sd(value))
    Col  <-  colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
    g <- ggplot(fields.long, aes(E_km, N_km, col=value)) +
      geom_point(size=1) +
      facet_grid(component+type~strata) +
      scale_colour_gradientn(colours = Col(15)) + theme_bw()
    ggsave(file.path(savedir, 'map_omegas.png'), g, width=9, height=6, units='in')
  }
  Report <- results$Report
  Index <- results$Index
  g <- ggplot(Index, aes(year, y=est, group=strata, fill=strata)) +
    geom_ribbon(aes(ymin=est-1.96*se, ymax=est+1.96*se), alpha=.5) +
    geom_line() + geom_point()+ theme_bw() +
    ylab('log abundance')
  ggsave(file.path(savedir, 'index.png'), g, width=7, height=5)
  ## Also create an index of the individual strata
  g <- ggplot(results$Index.strata, aes(year, y=est, group=strata, fill=strata)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
    geom_line() + geom_point()+ theme_bw() + ylab('log abundance')
  ggsave(file.path(savedir, 'index_strata.png'), g, width=7, height=5)
  Mapdetails <- make_map_info(Region, NN_Extrap=Spatial_List$NN_Extrap,
                              Extrapolation_List=Extrapolation_List)
  Mapdetails$Legend$x <- Mapdetails$Legend$x-70
  Mapdetails$Legend$y <- Mapdetails$Legend$y-45

  ## This was causing problems and not sure why. Will fix later.
  if(TmbData$n_c>1){
    Plot_factors(Report, results$ParHatList, Data=TmbData, SD=Opt$SD,
                 mapdetails_list=Mapdetails, plotdir=paste0(savedir, "/"))
  }
  Enc_prob <- plot_encounter_diagnostic(Report=Report,
                                        Data_Geostat=Data_Geostat,
                                        DirName=savedir)
  Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                               FileName_Phist="Posterior_Predictive-Histogram",
                               FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=savedir )

  MapDetails_List = make_map_info( "Region"=Region,
                                  "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
                                  "Extrapolation_List"=Extrapolation_List )
  MapDetails_List  <- Mapdetails
  ## Decide which years to plot
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'],
                 TmbData=TmbData, Report=Report, Q=Q, savedir=savedir,
                 MappingDetails=MapDetails_List[["MappingDetails"]],
                 PlotDF=MapDetails_List[["PlotDF"]],
                 MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
                 Xlim=MapDetails_List[["Xlim"]],
                 Ylim=MapDetails_List[["Ylim"]], FileName=savedir,
                 Year_Set=Year_Set, Years2Include=Years2Include,
                 Rotate=MapDetails_List[["Rotate"]],
                 Cex=MapDetails_List[["Cex"]],
                 Legend=MapDetails_List[["Legend"]],
                 zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  plot_anisotropy( FileName=paste0(savedir,"Aniso.png"), Report=Report,
                  TmbData=TmbData )
  ## Some built-in maps
  tmp <- c(1,2,3, 11, 12)
  if(results$Index$space[1]=='ST') tmp <- c(tmp, 6,7)
  Dens_xt = plot_maps(plot_set=tmp,
                      MappingDetails=MapDetails_List[["MappingDetails"]],
                      Report=Report, Sdreport=Opt$SD,
                      TmbData=TmbData,
                      PlotDF=MapDetails_List[["PlotDF"]],
                      MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
                      Xlim=MapDetails_List[["Xlim"]],
                      Ylim=MapDetails_List[["Ylim"]], FileName=paste0(savedir,'/'),
                      Year_Set=Year_Set, Years2Include=Years2Include,
                      Rotate=MapDetails_List[["Rotate"]],
                      Cex=MapDetails_List[["Cex"]],
                      Legend=MapDetails_List[["Legend"]],
                      zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
                      oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  Dens_DF = cbind( "Density"=as.vector(Dens_xt),
                  "Year"=Year_Set[col(Dens_xt)],
                  "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'],
                  "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

  Index = plot_biomass_index( DirName=savedir, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
  ##  pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )
  plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=savedir, Year_Set=Year_Set)

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
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/map_data_ratio_', ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=range(MatRatio, na.rm=TRUE),
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='log(Obs/Exp)', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
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
      PR2_i[i] <- (log(obs)-log(D_i[i]))/Report$SigmaM[gr]
    }
  }
  df <- cbind(Data_Geostat, PR1=PR1_i, PR2=PR2_i, positive=ifelse(Data_Geostat$Catch_KG>0,1,0))
  xlim <- range(df$Lon); ylim <- range(df$Lat)
  for(gr in 1:3){
    gt <- levels(Data_Geostat$Gear)[gr]
    g <- ggplot(subset(df, Gear==gt & positive==1), aes(Lon, Lat, size=abs(PR2), color=PR2>0))+
      geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
      scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
    ggsave(filename=paste0(savedir, '/Pearson_resid_catchrate_', gr, '.png'), plot=g,
           width=7, height=5)
    g <- ggplot(subset(df, Gear==gt & positive==0), aes(Lon, Lat, size=abs(PR1), color=PR1>0))+
      geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
      scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
    ggsave(filename=paste0(savedir, '/Pearson_resid_presence_', gr, '.png'), plot=g,
           width=7, height=5)
  }

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
## #### First get the Pearon resids for Trawl
## ## Sum across P1 for the first two strata then calculate R1 manually. I
## ## don't think there's another way to do this.
## exp_rate_xy <- 1-exp(-exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)))
## temp <- observed.catches.by.geartype('Trawl')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q1_xy[[1]] <- (obs_num_xy-exp_num_xy)/sqrt(exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS1
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('Acoustic_3-16')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q1_xy[[2]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS2
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('Acoustic_16-surface')
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
## #### First get the Pearon resids for Trawl
## ## Sum across P1 for the first two strata then calculate R1 manually. I
## ## don't think there's another way to do this.
## r1temp <- 1-exp(-exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)))
## ## For R2=exp(p1+p2)/r1
## exp_rate_xy <-
##   exp(apply(Report$P1_xcy[,-3,], c(1,3), sum)+apply(Report$P2_xcy[,-3,], c(1,3), sum))/r1temp
## temp <- observed.catches.by.geartype('Trawl')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q2_xy[[1]] <- (obs_num_xy-exp_num_xy)/sqrt(exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS1
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('Acoustic_3-16')
## exp_num_xy  <-  exp_rate_xy * temp$total_num_xy
## obs_num_xy  <-  temp$obs_rate_xy * temp$total_num_xy
## ## Now calculate Pearson residuals using binomial form:
## Q2_xy[[2]] <- (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(temp$total_num_xy-exp_num_xy)/temp$total_num_xy )
## #### Now ATS2
## exp_rate_xy <- 1-exp(-exp(Report$P1_xcy[,2,]))
## temp <- observed.catches.by.geartype('Acoustic_16-surface')
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





