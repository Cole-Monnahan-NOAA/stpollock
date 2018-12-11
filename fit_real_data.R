## File to run the fits to the real data

## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)


n_x <- 100 # number of knots
model <- 'combined'
space <- 'NS'
source("load_read_data.R")
Obj  <-  TmbList[["Obj"]]
Obj$env$beSilent()
Opt <- Optimize( obj=Obj, lower=TmbList[["Lower"]],
                  upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile,
                  bias.correct=FALSE, newtonsteps=1, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

## Manually calculate SE for the total biomass index. Since it's a sum of
## the three the derivatives are all 1 and so the SE is the sqrt of the sum
## of all of the variances and covariances. This feature is not coded
## into VAST yet so have to do it manually. also note that the order of
## the Index_cyl matrix in vector form is Index_11, Index_21, Index_31,
## Index_12,.. etc. This effects the subsetting below

est <- data.frame(par=names(Opt$SD$value), value=Opt$SD$value,
                  se=Opt$SD$sd)
tmp <- which(names(Opt$SD$value) %in% 'Index_cyl')
cov.index <- Opt$SD$cov[tmp,tmp]
index <- data.frame(year=years,
                    value=apply(Report$Index_cyl, 2, sum),
                    se=sqrt(sapply(1:nyr, function(i) {j=1:3+3*(i-1); sum(cov.index[j,j])})))
ggplot(index, aes(year, y=value)) +
  geom_ribbon(aes(ymin=value-1.96*se, ymax=value+1.96*se), fill=gray(.8)) +
                  geom_line() + theme_bw()


plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )
pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] )
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)
## Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
##                              FileName_Phist="Posterior_Predictive-Histogram",
##                              FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile )

MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
                                        # Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'],
               TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile,
               MappingDetails=MapDetails_List[["MappingDetails"]],
               PlotDF=MapDetails_List[["PlotDF"]],
               MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
               Xlim=MapDetails_List[["Xlim"]],
               Ylim=MapDetails_List[["Ylim"]], FileName=DateFile,
               Year_Set=Year_Set, Years2Include=Years2Include,
               Rotate=MapDetails_List[["Rotate"]],
               Cex=MapDetails_List[["Cex"]],
               Legend=MapDetails_List[["Legend"]],
               zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=1.8)

## plot_anisotropy( FileName=paste0(DateFile,"Aniso.png"), Report=Report,
##                 TmbData=TmbData )
## Dens_xt = plot_maps(plot_set=c(3),
##                 MappingDetails=MapDetails_List[["MappingDetails"]],
##                 Report=Report, Sdreport=Opt$SD,
##                 PlotDF=MapDetails_List[["PlotDF"]],
##                 MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
##                 Xlim=MapDetails_List[["Xlim"]],
##                 Ylim=MapDetails_List[["Ylim"]], FileName=DateFile,
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

Index = plot_biomass_index( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )

plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)
