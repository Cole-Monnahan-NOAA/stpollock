

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
Version <- "VAST_v5_3_0"

source("simulator.R")


process.results <- function(Opt, Obj, model, space, savedir){
  Report  <-  Obj$report()
  ParHat <- Opt$par #Obj$env$parList(Opt$par)
  Index <- calculate.index(Opt, Report, model, space)
  Save  <-  list(Index=Index, Opt=Opt, Report=Report, ParHat=ParHat, savedir=savedir)
  save(Save, file=paste0(savedir,"/Save.RData"))
  return(Save)
}
calculate.index <- function(Opt, Report, model, space){
  tmp <- which(names(Opt$SD$value) %in% 'Index_cyl')
  index <- data.frame(model=model, space=space,  year=years)
  if(model=='combined'){
  ## Manually calculate SE for the total biomass index. Since it's a sum of
  ## the three the derivatives are all 1 and so the SE is the sqrt of the sum
  ## of all of the variances and covariances. This feature is not coded
  ## into VAST yet so have to do it manually. also note that the order of
  ## the Index_cyl matrix in vector form is Index_11, Index_21, Index_31,
  ## Index_12,.. etc. This effects the subsetting below
    cov.index <- Opt$SD$cov[tmp,tmp]
    ## combined is all three strata summed
    index1 <- data.frame(index, strata='total',  est=apply(Report$Index_cyl[1:3,,], 2, sum),
                      se=sqrt(sapply(1:nyr, function(i) {j=1:3+3*(i-1);
                        sum(cov.index[j,j])})))
    ## sum the top two to get what the ATS sees
    index2 <- data.frame(index, strata='ats', est=apply(Report$Index_cyl[2:3,,], 2, sum),
                      se=sqrt(sapply(1:nyr, function(i) {j=(1:3+3*(i-1))[-1];
                        sum(cov.index[j,j])})))
    ## likewise the BTS is just the first strata
    index3 <- data.frame(index, strata='bts', est=Report$Index_cyl[1,,],
                      se=sqrt(sapply(1:nyr, function(i) {j=(1:3+3*(i-1))[1];
                        sum(cov.index[j,j])})))
    index <- rbind(index1,index2, index3)

  } else {
    ## chop off missing years for ATS case
    tmp2 <- which(min(years):max(years) %in% years)
    tmp <- tmp[tmp2]
    index <- data.frame(index, strata=model, est=Opt$SD$value[tmp], se=Opt$SD$sd[tmp])
  }
  index <- within(index, {lwr <- est-1.96*se; upr <- est+1.96*se})
  return(index)
}

plot.vastfit <- function(results){
  Report <- results$Report
  Index <- results$Index
  g <- ggplot(Index, aes(year, y=est, group=strata, fill=strata)) +
    geom_ribbon(aes(ymin=est-1.96*se, ymax=est+1.96*se), alpha=.5) +
    geom_line() + geom_point()+ theme_bw() + ylim(0, max(Index$est+2*Index$se))
  ggsave(file.path(savedir, 'index.png'), g, width=7, height=5)
  Mapdetails <- make_map_info(Region, NN_Extrap=Spatial_List$NN_Extrap,
                              Extrapolation_List=Extrapolation_List)
  Plot_factors(Report, results$ParHat, Data=TmbData, SD=Opt$SD,
               mapdetails_list=Mapdetails, plotdir=paste0(savedir, "/"))
  ParHat2 <- obj$env$parList()
  Enc_prob <- plot_encounter_diagnostic(Report=Report,
                                        Data_Geostat=Data_Geostat,
                                        DirName=savedir)
  Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                               FileName_Phist="Posterior_Predictive-Histogram",
                               FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=savedir )

  MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
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
                 oma=c(3.5,3.5,0,0), cex=1.8)

  ## plot_anisotropy( FileName=paste0(savedir,"Aniso.png"), Report=Report,
  ##                 TmbData=TmbData )
  ## Dens_xt = plot_maps(plot_set=c(3),
  ##                 MappingDetails=MapDetails_List[["MappingDetails"]],
  ##                 Report=Report, Sdreport=Opt$SD,
  ##                 PlotDF=MapDetails_List[["PlotDF"]],
  ##                 MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
  ##                 Xlim=MapDetails_List[["Xlim"]],
  ##                 Ylim=MapDetails_List[["Ylim"]], FileName=savedir,
  ##                 Year_Set=Year_Set, Years2Include=Years2Include,
  ##                 Rotate=MapDetails_List[["Rotate"]],
  ##                 Cex=MapDetails_List[["Cex"]],
  ##                 Legend=MapDetails_List[["Legend"]],
  ##                 zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
  ##                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  ## Dens_DF = cbind( "Density"=as.vector(Dens_xt),
  ##                 "Year"=Year_Set[col(Dens_xt)],
  ##                 "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'],
  ##                 "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

  Index = plot_biomass_index( DirName=savedir, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
  ##  pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )

  plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=savedir, Year_Set=Year_Set)
}
