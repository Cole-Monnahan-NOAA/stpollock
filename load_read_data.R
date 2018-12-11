Version <- "VAST_v5_3_0"

bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')
## The ats data is really high resolution so truncating this for now to
## make things faster and fix the mesh which is overly weighted to the ats
## data otherwise
ats <- ats[seq(1, nrow(ats), len=nrow(bts)),]

## bts <- subset(bts, year==2007)
## ats <- subset(ats, year==2007)

Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 50
n_x = 100
## Model settings
FieldConfig = c("Omega1"=3, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Delta1"=0, "Delta2"=0)
ObsModel = c(1,1)
Options =  c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=1, "Calculate_evenness"=0, "Calculate_effective_area"=1, "Calculate_Cov_SE"=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)
## Stratification for results
strata.limits <- data.frame('STRATA'="All_areas")
## Derived objects
Region = "Eastern_Bering_Sea"
## Save settings
DateFile = paste0(getwd(),'/VAST_output_real/')
dir.create(DateFile)
Record = list("Version"=Version,"Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
save( Record, file=file.path(DateFile,"Record.RData"))
capture.output( Record, file=paste0(DateFile,"Record.txt"))
TmbDir <- DateFile
DF_p1 = data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none')
DF_p2 = data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none')
DF_p3 = data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                   Vessel='none')
Data_Geostat = rbind( DF_p1, DF_p2, DF_p3 )
Extrapolation_List =
  make_extrapolation_info( Region=Region, strata.limits=strata.limits )
## Derived objects for spatio-temporal estimation
Spatial_List = make_spatial_info( grid_size_km=grid_size_km, n_x=n_x,
                                 Method=Method, Lon=Data_Geostat[,'Lon'],
                                 Lat=Data_Geostat[,'Lat'],
                                 Extrapolation_List=Extrapolation_List,
                                 DirPath=DateFile, Save_Results=FALSE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
c_iz = matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
              ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
# Add threshold
b_i = Data_Geostat[,'Catch_KG']
# Build data
TmbData = Data_Fn(Version=Version, FieldConfig=FieldConfig,
                  OverdispersionConfig=OverdispersionConfig,
                  RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                  b_i=b_i, a_i=Data_Geostat[,'AreaSwept_km2'],
                  v_i=as.numeric(Data_Geostat[,'Vessel'])-1,
                  s_i=Data_Geostat[,'knot_i']-1,
                  t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                  MeshList=Spatial_List$MeshList,
                  GridList=Spatial_List$GridList,
                  Method=Spatial_List$Method, Options=Options,
                  Aniso=FALSE)
Random = "generate"



