## Central file to fit the pollock data using the modified VAST model

library(devtools)
## remove.packages('VAST')
## install_github('James-Thorson/VAST')
library(VAST)
library(TMB)
library(maps)
library(mapdata)

## I stripped these out of the VAST single species example
Version = get_latest_version( package="VAST" )
## Spatial settings
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
n_x = 100   # Specify number of stations (a.k.a. "knots")
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
## FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,0)
DateFile = paste0(getwd(),'/VAST_output/')
dir.create(DateFile)

Kmeans_Config <- list('randomseed'=1, 'nstart'=1, 'iter.max'=200)

strata.limits <- data.frame('STRATA'="All_areas")
Region <- "Eastern_Bering_Sea"
data( EBS_pollock_data, package="FishStatsUtils" )
Data_Geostat = data.frame( "Catch_KG"=EBS_pollock_data[,'catch'],
                          "Year"=EBS_pollock_data[,'year'],
                          "Vessel"="missing", "AreaSwept_km2"=0.01,
                          "Lat"=EBS_pollock_data[,'lat'],
                          "Lon"=EBS_pollock_data[,'long'], "Pass"=0)
Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )
Spatial_List = make_spatial_info( grid_size_km=grid_size_km, n_x=n_x,
                                 Method=Method, Lon=Data_Geostat[,'Lon'],
                                 Lat=Data_Geostat[,'Lat'],
                                 Extrapolation_List=Extrapolation_List,
                                 randomseed=Kmeans_Config[["randomseed"]],
                                 nstart=Kmeans_Config[["nstart"]],
                                 iter.max=Kmeans_Config[["iter.max"]],
                                 DirPath=DateFile, Save_Results=FALSE )
# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
Options =  c("SD_site_density"=0, "SD_site_logdensity"=0,
             "Calculate_Range"=1, "Calculate_evenness"=0,
             "Calculate_effective_area"=1, "Calculate_Cov_SE"=0,
             'Calculate_Synchrony'=0, 'Calculate_Coherence'=0,
             'normalize_GMRF_in_CPP'=TRUE)

TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig, "ObsModel"=ObsModel,
                  "c_i"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method, "Options"=Options )
TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile,
                  "Version"=Version, "RhoConfig"=RhoConfig,
                  "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=1, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )
pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] )
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)
Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                             FileName_Phist="Posterior_Predictive-Histogram",
                             FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile )

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

plot_anisotropy( FileName=paste0(DateFile,"Aniso.png"), Report=Report,
                TmbData=TmbData )
Dens_xt = plot_maps(plot_set=c(3),
                MappingDetails=MapDetails_List[["MappingDetails"]],
                Report=Report, Sdreport=Opt$SD,
                PlotDF=MapDetails_List[["PlotDF"]],
                MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
                Xlim=MapDetails_List[["Xlim"]],
                Ylim=MapDetails_List[["Ylim"]], FileName=DateFile,
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

Index = plot_biomass_index( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )
plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)
