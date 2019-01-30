### This file is meant to be sourced given some global options, resulting
### in inputs ready for use in VAST. There are two main options: model type
### (ATS only, BTS only, or combined) and then three versions of spatial
### complexity (no space [NS], space [S], spatiotemporal [ST]). These
### values trigger different configurations in the code below
stopifnot(model %in% c('ats', 'bts', 'combined'))
stopifnot(space %in% c('NS', 'S', 'ST'))

## Default to suppress messages to cleanup output
if(!exists('silent')) silent <- TRUE
silent.fn <- function(expr){
  if(silent) suppressMessages(expr) else expr
}

### Step 1: Load in the real data
DF1 <- data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none', depth=bts$depth)
DF2 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth)
DF3 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth)

### Step 2: Configure the spatial factors which depend on inputs
n_f <- ifelse(model=='combined', 3,1) # number of factors to use
## This puts a FA on beta1 and beta2, which means I need to set Rho
## accordingly below
FieldConfig <- matrix(c("Omega1"=ifelse(space=='NS', 0,n_f),
                        "Epsilon1"=ifelse(space=='ST', n_f,0),
                        "Beta1"=n_f,
                        "Omega2"=0,
                        "Epsilon2"=0, "Beta2"=n_f), ncol=2 )
### Rho config= 0: each year as fixed effect; 1: each year as random
### following IID distribution; 2: each year as random following a random
### walk; 3: constant among years as fixed effect; 4: each year as random
### following AR1 process
## For now using IID for combined model and temporal on ATS/BTS since
## missing years there.
x <- switch(model, combined=1, ats=4, bts=4)
RhoConfig <- c("Beta1"=x, "Beta2"=x, "Epsilon1"=0, "Epsilon2"=0)


### Step 3: Setup VAST inputs which are constant for the models
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 50
## Model settings
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
ObsModel <- c(1,1)
Options <-  c("SD_site_density"=0, "SD_site_logdensity"=0,
              "Calculate_Range"=1, "Calculate_evenness"=0,
              "Calculate_effective_area"=1, "Calculate_Cov_SE"=0,
              'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)
## Stratification for results
strata.limits <- data.frame('STRATA'="All_areas")
Region <- "Eastern_Bering_Sea"
silent.fn(Extrapolation_List <-
            make_extrapolation_info(Region=Region, strata.limits=strata.limits))

## Derived objects
## Save settings
## savedir <- paste0(getwd(),'/VAST_output_real/')
dir.create(savedir, showWarnings=FALSE)
## ## Copy over the DLL so I don't have to compile each time.
trash <- file.copy(file.path('models', paste0(Version, '.cpp')),
          to=file.path(savedir, paste0(Version, '.cpp')))
trash <- file.copy(file.path('models', paste0(Version, '.dll')),
          to=file.path(savedir, paste0(Version, '.dll')))
trash <- file.copy(file.path('models', paste0(Version, '.o')),
          to=file.path(savedir, paste0(Version, '.o')))
Record <- list("Version"=Version,"Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
save( Record, file=file.path(savedir,"Record.RData"))
capture.output( Record, file=paste0(savedir,"/Record.txt"))

### Step 4: Construct VAST model based on inputs and data
Q_ik <- NULL ## catchability covariates, updated below for combined model?
if(model=='combined'){
  Data_Geostat <- rbind( DF1, DF2, DF3 )
  c_iz <- matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
                 ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
  ## Q_ik <- cbind(ifelse(Data_Geostat$Gear=='Trawl', 1, 0),
  ##               ifelse(Data_Geostat$Gear=='Trawl', 0, 1))
} else if(model=='ats'){
  ## For this one sum across the two strata to create a single one, akin to
  ## what they'd do without the BTS
  Data_Geostat <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2+ats$strata3, depth=ats$depth,
                   Gear='Acoustic_3-surface', AreaSwept_km2=1,
                   Vessel='none')
  c_iz <- rep(0, nrow(Data_Geostat))
  years <- sort(unique(ats$year))
} else if(model=='bts'){
  Data_Geostat <- DF1
  years <- sort(unique(bts$year))
  c_iz <- rep(0, nrow(Data_Geostat))
}
years <- sort(unique(Data_Geostat$Year))
nyr <- length(years)
## Derived objects for spatio-temporal estimation
silent.fn(Spatial_List <-
     make_spatial_info(grid_size_km=grid_size_km, n_x=n_x, Method=Method,
                       Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],
                       Extrapolation_List=Extrapolation_List, DirPath=savedir,
                       Save_Results=FALSE ))
Data_Geostat <- cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
silent.fn(XX <- (FishStatsUtils::format_covariates(
    Lat_e = Data_Geostat$Lat,
    Lon_e = Data_Geostat$Lon,
    t_e = Data_Geostat$Year,
    Cov_ep = Data_Geostat[,'depth'],
    Extrapolation_List = Extrapolation_List,
    Spatial_List = Spatial_List, FUN = mean,
    na.omit = "time-average")))
##dimnames(X_xtp)[[1]] <- dimnames(covsperknot$Cov_xtp)[[1]]
X_xtp <- XX$Cov_xtp

## Build data and object for first time
TmbData <- Data_Fn(Version=Version, FieldConfig=FieldConfig,
                  OverdispersionConfig=OverdispersionConfig,
                  RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                  b_i=Data_Geostat[,'Catch_KG'],
                  a_i=Data_Geostat[,'AreaSwept_km2'],
                  ## v_i=as.numeric(Data_Geostat[,'Vessel'])-1,
                  s_i=Data_Geostat[,'knot_i']-1,
                  t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                  MeshList=Spatial_List$MeshList,
                  GridList=Spatial_List$GridList,
                  Q_ik=Q_ik,
                  X_xtp=X_xtp,
                  Method=Spatial_List$Method, Options=Options,
                  Aniso=FALSE)
TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=savedir,
                         Version=Version,  RhoConfig=RhoConfig,
                         loc_x=Spatial_List$loc_x, Method=Method,
                         TmbDir='models', Random="generate")
## Tweak the Map based on inputs
Map <- TmbList0$Map
Params <- TmbList0$Parameters
if(model=='combined'){
  ## Assume that the two ATS strata have the same observation error
  Map$logSigmaM <- factor( cbind( c(1,2,2), NA, NA) )
  ##  Map$beta1_ct <- factor(rep(1, 30))
  ## Carefully build the catchability
  ## Params$lambda1_k <- c(0,0)
  ## ## Leave the first fixed otherwise confounded with the betas (RIGHT?)
  ## Map$lambda1_k <- as.factor(c(NA,1))
  ## Map$lambda2_k <- as.factor(c(NA,NA))
} else if(model=='ats'){
## Estimate a single parameter for the second LP regardless of model. Need
## to be careful to not estimate years without data in the 'ats' case where
## VAST already uses a map with NA for missing years.
  Map$beta2_ct[which(!is.na(Map$beta2_ct))] <- 1
  Map$beta2_ct <- droplevels(as.factor(Map$beta2_ct))
} else if(model=='bts'){
  ## This has no NA b/c all years represented in the data
  Map$beta2_ct <- factor(rep(1, length(Params$beta2_ct)))
}
## Turn off the depth effect for the second LP and initialize at 0
Params$gamma1_ctp <- Params$gamma2_ctp <- Params$gamma1_ctp*0
## Map$gamma1_ctp  <- factor(rep(1,length(Params$gamma1_ctp)))
Map$gamma2_ctp  <- factor(Params$gamma2_ctp*NA)

## if(space=='ST' & model =='combined'){
##   ## turn off estimation of factor analysis and just do diagonal (for now)
##   n_f <- 3; tmp <- diag(1:n_f, nrow=3, ncol=n_f)
##   lvec <- tmp[lower.tri(tmp, TRUE)] # init values
##   Map$L_epsilon1_z <- factor(ifelse(lvec==0, NA, lvec))
##   Params$L_epsilon_z <- lvec
## }

## Rebuild with the new mapping stuff
TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=savedir,
                        Version=Version,  RhoConfig=RhoConfig,
                        loc_x=Spatial_List$loc_x, Method=Method,
                        Param=Params, TmbDir='models', Random='generate',
                         Map=Map)
Obj  <-  TmbList[["Obj"]]
Obj$env$beSilent()

## bundle together some of the inputs that will be needed later for
## plotting and such that aren't included in the standard VAST output
loc <- data.frame(Spatial_List$MeshList$isotropic_mesh$loc[,-3])
names(loc) <- c('E_km', 'N_km')
Inputs <- list(loc=loc)


silent.fn(plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List,
          Data_Geostat=Data_Geostat, PlotDir=paste0(savedir,"/") ))

