## source this to get the data ready for the models in compare_models.R
## script


library(TMB)
library(TMBhelper)
library(INLA)
library(reshape2)
library(ggplot2)
library(VAST)


layers <- as.matrix(read.csv('data/layers.csv', header=FALSE, sep=','))
hauls <- read.csv('data/hauls.csv')
SA <- layers; dimnames(SA) <- NULL
## subset down to a single year to keep it simple
layers <- layers[which(hauls$year==2006),]
hauls <- hauls[which(hauls$year==2006),]
ntows <- nrow(hauls)
## From Stan: first two columns are the ADZ so skip them. 3rd is 0.5-0.75m
## and should assume this is 'h2' or 'h' in the paper. More specifically:
## the first 20 columns are vertical layers from 0 to 5m, every 0.25m, rest
## of the columns are in 1m intervals
layer.widths <- c(seq(.25,5, by=.25), 6+1:(ncol(SA)-20))
atstar <- at1 <- at2<- rep(NA, ntows)
EFH <- 14; iEFH <- which(layer.widths==EFH)
for (i in 1:ntows){
  atstar[i] <- SA[i, 3] # layer just above the ADZ (was sum_SA2)
  at1[i] <- sum(SA[i, 3:iEFH]) # ADZ to EFH (was sum_SA1)
  at2[i] <- sum(SA[i, (iEFH+1):ncol(SA)]) # EFH to surface
}
bt <- hauls$pred_sa
## The three strata data
Y <- cbind(log(bt), log(at1), log(at2))


## prepare the factor analysis component
tmp <- diag(1:n_f, nrow=3, ncol=n_f)
lvec <- tmp[lower.tri(tmp, TRUE)] # init values
## This map will turn off estimation of off diagonals but estimate the
## diagonal which turns off the FA part of the model
fa.off <- factor(ifelse(lvec==0, NA, lvec))

## prepare the geospatial meshxb
mesh <-  inla.mesh.create( hauls[,c('s_long', 's_lat')])
spde <- inla.spde2.matern( mesh )


### The VAST inputs for the same data
### Fit VAST to this same data subset
Version = "VAST_v4_0_0"
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
DateFile = paste0(getwd(),'/VAST_output/')
dir.create(DateFile)
Record = list("Version"=Version,"Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
save( Record, file=file.path(DateFile,"Record.RData"))
capture.output( Record, file=paste0(DateFile,"Record.txt"))
TmbDir <- DateFile
DF_p1 = data.frame( Lat=hauls$s_lat, Lon=hauls$s_long, Year=2006,
                   Catch_KG=bt, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none')
DF_p2 = data.frame( Lat=hauls$s_lat, Lon=hauls$s_long, Year=2006,
                   Catch_KG=at1, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none')
DF_p3 = data.frame( Lat=hauls$s_lat, Lon=hauls$s_long, Year=2006,
                   Catch_KG=at2, Gear='Acoustic_16-surface', AreaSwept_km2=1,
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
# Build and run model
## Build model
#source( "C:/Users/James.Thorson/Desktop/Project_git/VAST/R/Data_Fn.R" )
#c_iz = as.numeric(Data_Geostat[,'Gear']) - 1
c_iz = matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
              ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
# Add threshold
b_i = Data_Geostat[,'Catch_KG']
## b_i = ifelse( b_i<100, 0, b_i )
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
                  Method=Spatial_List$Method, Options=Options )
Random = "generate"
