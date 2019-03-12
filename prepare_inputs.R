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
source("load_data.R")

### Step 2: Configure the spatial factors which depend on inputs
n_f <- ifelse(model=='combined', 3,1) # number of factors to use
n_eps <- 2
## This puts a FA on beta1 and beta2, which means I need to set Rho
## accordingly below
FieldConfig <- matrix(c("Omega1"=ifelse(space=='NS', 0,n_f),
                        "Epsilon1"=ifelse(space=='ST',n_eps,0),
                        "Beta1"="IID",
                        "Omega2"=ifelse(space=='NS', 0, n_f),
                        "Epsilon2"=ifelse(space=='ST', 0,0),
                        "Beta2"="IID"), ncol=2 )
### Rho config= 0: each year as fixed effect; 1: each year as random
### following IID distribution; 2: each year as random following a random
### walk; 3: constant among years as fixed effect; 4: each year as random
### following AR1 process
## For now using IID for combined model and temporal on ATS/BTS since
## missing years there.
x <- switch(model, combined=4, ats=4, bts=4)
RhoConfig <- c("Beta1"=x, "Beta2"=3, "Epsilon1"=ifelse(space=='ST', 4,0), "Epsilon2"=0)


### Step 3: Setup VAST inputs which are constant for the models
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 50
## Model settings
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
ObsModel <- c(1,1)
Options <-  c("SD_site_density"=0, "SD_site_logdensity"=1,
              "Calculate_Range"=1, "Calculate_evenness"=0,
              "Calculate_effective_area"=1, "Calculate_Cov_SE"=1,
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
  ## Data_Geostat <- subset(Data_Geostat, Year < 2011)
  ##  tmp <- Data_Geostat
  ##  tmp$Year <- tmp$Year-4
  ##  Data_Geostat <- rbind(Data_Geostat, tmp)
  ##  tmp$Year <- tmp$Year-4
  ##  Data_Geostat <- rbind(Data_Geostat, tmp)
  c_iz <- matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
                 ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
  ## c_iz[,2] <- NA
  Q_ik <- matrix(ifelse(Data_Geostat$Gear=='Trawl', 1, 0), ncol=1)
} else if(model=='ats'){
  ## For this one sum across the two strata to create a single one, akin to
  ## what they'd do without the BTS
  Data_Geostat <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                             Catch_KG=ats$strata2+ats$strata3, depth=ats$depth,
                             depth2=ats$depth2,
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

### Derived objects for spatio-temporal estimation
silent.fn(Spatial_List  <-
  make_spatial_info(grid_size_km=grid_size_km, n_x=n_x,
                    Method=Method, Lon=Data_Geostat[,'Lon'],
                    Lat=Data_Geostat[,'Lat'],
                    ## According to Jim this will make the grid uniform with respect to the
                    ## extrapolation region. This helps avoid the grid being driven by the ATS
                    ## which as more points than the BTS. But it breaks
                    ## plotting code so turned off for now.
                    ## LON_intensity=Extrapolation_List$Data_Extrap[which(Extrapolation_List$Data_Extrap$Include==1),'Lon'],
                    ## LAT_intensity=Extrapolation_List$Data_Extrap[which(Extrapolation_List$Data_Extrap$Include==1),'Lat'],
                    Extrapolation_List=Extrapolation_List,
                    DirPath=savedir, Save_Results=FALSE ))
## silent.fn(Spatial_List <-
##             make_spatial_info(grid_size_km=grid_size_km, n_x=n_x, Method=Method,
##                               Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],
##                               Extrapolation_List=Extrapolation_List, DirPath=savedir,
##                               Save_Results=FALSE ))
Data_Geostat <- cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
silent.fn(XX <- (FishStatsUtils::format_covariates(
                                   Lat_e = Data_Geostat$Lat,
                                   Lon_e = Data_Geostat$Lon,
                                   t_e = Data_Geostat$Year,
                                   Cov_ep = Data_Geostat[,c('depth', 'depth2')],
                                   Extrapolation_List = Extrapolation_List,
                                   Spatial_List = Spatial_List, FUN = mean,
                                   na.omit = "time-average")))

## Build data and object for first time
TmbData <- Data_Fn(Version=Version, FieldConfig=FieldConfig,
                   OverdispersionConfig=OverdispersionConfig,
                   RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                   b_i=Data_Geostat[,'Catch_KG'],
                   a_i=Data_Geostat[,'AreaSwept_km2'],
                   ## v_i=as.numeric(Data_Geostat[,'Vessel'])-1,
                   v_i=1:nrow(Data_Geostat),
                   s_i=Data_Geostat[,'knot_i']-1,
                   t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                   MeshList=Spatial_List$MeshList,
                   GridList=Spatial_List$GridList,
                   Q_ik=Q_ik,
                   X_xtp=NULL,#XX$Cov_xtp,
                   Method=Spatial_List$Method, Options=Options,
                   Aniso=FALSE)
TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=savedir,
                         Version=Version,  RhoConfig=RhoConfig,
                         loc_x=Spatial_List$loc_x, Method=Method,
                         TmbDir='models', Random="generate")
## Tweak the Map based on inputs
Map <- TmbList0$Map
Params <- TmbList0$Parameters
Params$Beta_mean2_c <- Params$Beta_mean2_c+5
if(model=='combined'){
  Params$L_beta1_z <- c(.2,.3,.5)
  Params$L_beta2_z <- c(.6,.3,1)
  Params$logSigmaM[1:3] <- c(.6,.7,.8)
  ## Map$lambda1_k <- Map$lambda2_k <- factor(NA)
} else {
  Params$L_beta1_z <- .4
  Params$L_beta2_z <- .4
}
Params$logkappa1 <- Params$logkappa2 <- -5
if(space=='ST'){
  if(length(Params$Beta_rho1_f)!=3) stop('problem with beta_rho1')
  Map$Beta_rho1_f <- factor(c(2,2,2))
}

## Params$beta2_ft <- Params$beta2_ft+5
## if(model=='combined'){
##   ## Assume that the two ATS strata have the same observation error
##  ## Map$logSigmaM <- factor( cbind( c(1,2,2), NA, NA) )
##   ##  Map$beta1_ct <- factor(rep(1, 30))
##   ## Carefully build the catchability
##   ## Params$lambda1_k <- c(0,0)
##   ## ## Leave the first fixed otherwise confounded with the betas (RIGHT?)
##   ## Map$lambda1_k <- as.factor(c(NA,1))
##   ## Map$lambda2_k <- as.factor(c(NA,NA))
## } else if(model=='ats'){
## ## Estimate a single parameter for the second LP regardless of model. Need
## ## to be careful to not estimate years without data in the 'ats' case where
## ## VAST already uses a map with NA for missing years.
##   Map$beta2_ct[which(!is.na(Map$beta2_ct))] <- 1
##   Map$beta2_ct <- droplevels(as.factor(Map$beta2_ct))
## } else if(model=='bts'){
##   ## This has no NA b/c all years represented in the data
##   Map$beta2_ct <- factor(rep(1, length(Params$beta2_ct)))
## }
## Set depth and depth2 coefficients to be constant across years and strata
## but affecting p1 and p2
## Params$gamma1_ctp <- Params$gamma2_ctp <- Params$gamma1_ctp*0
## tmp <- Params$gamma1_ctp
## tmp[,,1] <- 1; tmp[,,2] <- 2 # depth and depth2 are separate
## Map$gamma1_ctp <- Map$gamma2_ctp <- factor(tmp)

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
                        Param=Params, TmbDir='models',
                        Random='generate',
                        ##Random=c('beta1_ft'),
                        Map=Map)
Obj  <-  TmbList[["Obj"]]
Obj$env$beSilent()

## bundle together some of the inputs that will be needed later for
## plotting and such that aren't included in the standard VAST output
loc <- data.frame(Spatial_List$MeshList$isotropic_mesh$loc[,-3])
names(loc) <- c('E_km', 'N_km')
Inputs <- list(loc=loc, loc_x=data.frame(knot_x=1:n_x, Spatial_List$loc_x))


silent.fn(plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List,
                    Data_Geostat=Data_Geostat, PlotDir=paste0(savedir,"/") ))


## Some custom maps of the data properties
## Plot log average catch in grid
mdl <- make_map_info( "Region"=Region,
                     "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
                     "Extrapolation_List"=Extrapolation_List )
mdl$Legend$x <- mdl$Legend$x-70
mdl$Legend$y <- mdl$Legend$y-45
Year_Set <- sort(unique(Data_Geostat$Year))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
MatDat <- log(tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                     FUN=mean, na.rm=TRUE))
## Some grids have only zero observations
MatDat[is.infinite(MatDat)]  <-  NA
## Use consistent zlim for all three data types
zlim <- range(MatDat, na.rm=TRUE)
if(model=='combined'){
  message('Making data maps by gear type...')
  for(ii in 1:3){
    PlotMap_Fn(MappingDetails=mdl$MappingDetails,
               Mat=MatDat[,ii,Years2Include,drop=TRUE],
               PlotDF=mdl$PlotDF,
               MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_data_avg_', ii),
               Year_Set=Year_Set[Years2Include],
               Legend=mdl$Legend, zlim=zlim,
               mfrow = c(ceiling(sqrt(length(Years2Include))),
                         ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
               textmargin='Log avg catches', zone=mdl$Zone, mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
  }

  ## Plot percentage 0's
  MatDat <- tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                   FUN=function(x) mean(x>0, na.rm=TRUE))
  for(ii in 1:3){
    PlotMap_Fn(MappingDetails=mdl$MappingDetails,
               Mat=MatDat[,ii,Years2Include,drop=TRUE],
               PlotDF=mdl$PlotDF,
               MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
               FileName=paste0(savedir, '/map_data_pres_', ii),
               Year_Set=Year_Set[Years2Include],
               Legend=mdl$Legend, zlim=c(0,1),
               mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
               textmargin='Presence', zone=mdl$Zone, mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
  }


  ## Plot log average catch in BTS divided by ATS 3-16. This shouldn't be
  ## possible b/c the BTS also includes the ATS data. Although this ignores
  ## catchability.
  MatDat <- (tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                    FUN=mean, na.rm=TRUE))
  MatDat <- (MatDat[,2,]/MatDat[,1,])>1
  ## Some grids have only zero observations
  MatDat[is.infinite(MatDat) | MatDat==0]  <-  NA
  PlotMap_Fn(MappingDetails=mdl$MappingDetails,
             Mat=MatDat[,Years2Include],
             PlotDF=mdl$PlotDF, zlim=c(0,1),
             MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
             FileName=paste0(savedir, '/map_data_ratio'),
             Year_Set=Year_Set[Years2Include],
             Legend=mdl$Legend,
             mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
             textmargin='Ratio log(ATS)/log(BTS)) avg catches', zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
             oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
}

