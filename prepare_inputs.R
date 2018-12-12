## This file is meant to be sourced given some global options, resulting in
## inputs ready for use in VAST.

## number of factors for ST
if(space=='ST'){
  ## was having convergence issues with ST model so try a simpler ST
  ## component by turning off estimation of L_espilon1_z below.
  neps <- ifelse(model=='combined', 3, 1)
} else {
  neps <- 0
}

## Factors for space
n_f <- ifelse(model=='combined', 3,1)
FieldConfig <- c("Omega1"=n_f, "Epsilon1"=neps, "Omega2"=0, "Epsilon2"=0)



## Load in the real data
bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')
## The ats data is really high resolution so truncating this for now to
## make things faster and fix the mesh which is overly weighted to the ats
## data otherwise
ats <- ats[seq(1, nrow(ats), len=nrow(bts)),]

## Setup VAST inputs
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 50

## Model settings
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
ObsModel <- c(1,1)
Options <-  c("SD_site_density"=0, "SD_site_logdensity"=0,
              "Calculate_Range"=1, "Calculate_evenness"=0,
              "Calculate_effective_area"=1, "Calculate_Cov_SE"=0,
              'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)
## Stratification for results
strata.limits <- data.frame('STRATA'="All_areas")
## Derived objects
Region <- "Eastern_Bering_Sea"
## Save settings
savedir <- paste0(getwd(), '/VAST_output_', model, "_", space)
## savedir <- paste0(getwd(),'/VAST_output_real/')
dir.create(savedir, showWarnings=FALSE)
## Copy over the DLL so I don't have to compile each time.
trash <- file.copy(file.path('models', paste0(Version, '.cpp')),
          to=file.path(savedir, paste0(Version, '.cpp')))
trash <- file.copy(file.path('models', paste0(Version, '.dll')),
          to=file.path(savedir, paste0(Version, '.dll')))
trash <- file.copy(file.path('models', paste0(Version, '.o')),
          to=file.path(savedir, paste0(Version, '.o')))

Record <- list("Version"=Version,"Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
save( Record, file=file.path(savedir,"Record.RData"))
capture.output( Record, file=paste0(savedir,"/Record.txt"))
TmbDir <- savedir
DF_p1 <- data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none')
DF_p2 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none')
DF_p3 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                   Vessel='none')

if(model=='combined'){
  Data_Geostat <- rbind( DF_p1, DF_p2, DF_p3 )
  c_iz <- matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
              ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
} else if(model=='ats'){
  ## For this one sum across the two strata to create a single one, akin to
  ## what they'd do without the BTS
  Data_Geostat <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2+ats$strata3,
                   Gear='Acoustic_3-surface', AreaSwept_km2=1,
                   Vessel='none')
  c_iz <- rep(0, nrow(Data_Geostat))
  years <- sort(unique(ats$year))
  nyr <- length(years)
} else if(model=='bts'){
  Data_Geostat <- DF_p1
  years <- sort(unique(bts$year))
  nyr <- length(years)
  c_iz <- rep(0, nrow(Data_Geostat))
} else {
  stop("invalid model type")
}

Extrapolation_List =
  make_extrapolation_info( Region=Region, strata.limits=strata.limits )
## Derived objects for spatio-temporal estimation
Spatial_List <- make_spatial_info( grid_size_km=grid_size_km, n_x=n_x,
                                 Method=Method, Lon=Data_Geostat[,'Lon'],
                                 Lat=Data_Geostat[,'Lat'],
                                 Extrapolation_List=Extrapolation_List,
                                 DirPath=savedir, Save_Results=FALSE )
Data_Geostat <- cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Add threshold
b_i <- Data_Geostat[,'Catch_KG']
# Build data
TmbData <- Data_Fn(Version=Version, FieldConfig=FieldConfig,
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
Random <- "generate"
TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=savedir,
                         Version=Version,  RhoConfig=RhoConfig,
                         loc_x=Spatial_List$loc_x, Method=Method,
                         TmbDir='models', Random="generate")


Map <- TmbList0$Map
Params <- TmbList0$Parameters
## Fix SigmaM for all surveys to be equal
if(model=='combined'){
  Map$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
  ##  Map$beta1_ct <- factor(rep(1, 30))
}
## Estimate a single parameter for the second LP regardless of model
Map$beta2_ct <- factor(rep(1, length(Params$beta2_ct)))

if(space == 'NS'){
  ## turn off estimation of space
  Map$logkappa1 <- factor(NA); Params$logkappa1 <- 5
}
if(space=='ST' & model =='combined'){
  ## turn off estimation of factor analysis and just do diagonal (for now)
  n_f <- 3; tmp <- diag(1:n_f, nrow=3, ncol=n_f)
  lvec <- tmp[lower.tri(tmp, TRUE)] # init values
  Map$L_epsilon1_z <- factor(ifelse(lvec==0, NA, lvec))
  Params$L_epsilon_z <- lvec
}
TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=savedir,
                        Version=Version,  RhoConfig=RhoConfig,
                        loc_x=Spatial_List$loc_x, Method=Method,
                        Param=Params, TmbDir=TmbDir, Random='generate',
                         Map=Map)

Obj  <-  TmbList[["Obj"]]
Obj$env$beSilent()

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List,
          Data_Geostat=Data_Geostat, PlotDir=paste0(savedir,"/") )

