### This file is meant to be sourced given some global options, resulting
### in inputs ready for use in VAST. There are two main options: model type
### (ATS only, BTS only, or combined) and then three versions of spatial
### complexity (no space [NS], space [S], spatiotemporal [ST]). These
### values trigger different configurations in the code below


## Setup default configuration if not specified in input control list
finescale <- ifelse(is.null(control$finescale), FALSE, control$finescale)
## default is to estimate both lambdas
fixlambda <- ifelse(is.null(control$fixlambda), 0, control$fixlambda)
filterdata <- ifelse(is.null(control$filterdata), TRUE, control$filterdata)
filteryears <- ifelse(is.null(control$filteryears), FALSE, control$filteryears)
simdata <- ifelse(is.null(control$simdata), FALSE, control$simdata)
combinedoff <- ifelse(is.null(control$combinedoff), FALSE, control$combinedoff)
make_plots <- ifelse(is.null(control$make_plots), FALSE, control$make_plots)
silent.console <- ifelse(is.null(control$silent.console), TRUE, control$silent.console)
silent  <- ifelse(is.null(control$silent ), TRUE, control$silent )
temporal <- ifelse(is.null(control$temporal), 2, control$temporal)
beta1temporal <- ifelse(is.null(control$beta1temporal), TRUE, control$beta1temporal)
beta2temporal <- ifelse(is.null(control$beta2temporal), TRUE, control$beta2temporal)
kappaoff <- ifelse(is.null(control$kappaoff), 12, control$kappaoff)
seed <- ifelse(is.null(control$seed), 9999, control$seed)
n_x <- ifelse(is.null(control$n_x), 50, control$n_x)
model <- ifelse(is.null(control$model), 'combined', control$model)

set.seed(seed)

## These depend on the model and spatial setup
if(model != 'combined'){
  n_omega1 <- ifelse(is.null(control$n_omega1), 1, control$n_omega1)
  n_omega2 <- ifelse(is.null(control$n_omega2), 1, control$n_omega2)
  n_eps1 <- ifelse(is.null(control$n_eps1), 1, control$n_eps1)
  n_eps2 <- ifelse(is.null(control$n_eps2), 1, control$n_eps2)
} else {
  n_omega1 <- ifelse(is.null(control$n_omega1), 2, control$n_omega1)
  n_omega2 <- ifelse(is.null(control$n_omega2), 0, control$n_omega2)
  n_eps1 <- ifelse(is.null(control$n_eps1), 2, control$n_eps1)
  n_eps2 <- ifelse(is.null(control$n_eps2), 0, control$n_eps2)
}

## if(space!='ST'){
##   ## spatial only
##   n_eps1 <- n_eps2 <- 0
##   if(space=='NS'){
##     n_omega1 <- n_omega2 <- 0
##   }
## }

#stopifnot(temporal %in% c(2,4)) ## RW or AR1 only
stopifnot(model %in% c('ats', 'bts', 'combined'))
## stopifnot(space %in% c('NS', 'S', 'ST'))
## Be careful since these n's can be equal to "IID" which does weird things
## with logical statements
if(any(n_eps1>0, n_eps2>0)){
  space <- 'ST'
} else if(any(n_omega1>0,n_omega2>0)){
  space <- 'S'
} else {
  space <- 'NS'
}



## Default to suppress messages to cleanup output
silent.fn <- function(expr){
  if(silent.console) suppressMessages(expr) else expr
}

### Step 1: Load in the real data if not doing simulation
if(!simdata){
  source("load_data.R")
} else {
  message("Using simulated data DF1, DF2, DF3 in global workspace")
}

### Step 2: Configure the spatial factors which depend on inputs
FieldConfig <- matrix(c("Omega1"= n_omega1,
                        "Epsilon1"=n_eps1,
                        "Beta1"="IID",
                        "Omega2"=n_omega2,
                        "Epsilon2"=n_eps2,
                        "Beta2"="IID"), ncol=2 )
### Rho config= 0: each year as fixed effect; 1: each year as random
### following IID distribution; 2: each year as random following a random
### walk; 3: constant among years as fixed effect; 4: each year as random
### following AR1 process
## For now using IID for combined model and temporal on ATS/BTS since
## missing years there.
RhoConfig <- c("Beta1"=ifelse(beta1temporal, temporal, 3),
               "Beta2"=ifelse(beta2temporal, temporal, 3),
               "Epsilon1"=ifelse(n_eps1>0, temporal, 0),
               "Epsilon2"=ifelse(n_eps2>0, temporal, 0))

### Step 3: Setup VAST inputs which are constant for the models
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 50
## Model settings
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
ObsModel <- c(1,1)
Options <-  c("SD_site_density"=0, "SD_site_logdensity"=1,
              "Calculate_Range"=1, "Calculate_evenness"=0,
              "Calculate_effective_area"=1, "Calculate_Cov_SE"=1,
              'Calculate_Synchrony'=0, 'Calculate_Coherence'=0,
              'Calculate_proportion'=1)
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
  ## This is a switch to turn off the combined part and revert back to
  ## standard multivariate model. For testing only.
  if(combinedoff){ c_iz[,2] <- NA; warning('turned off combined part')}
  if(fixlambda== -1){
    message('turning on annual catchability')
    ## Annual coefficient, create model matrix with sum to zero contrasts
    ## to avoid confounding with betas
    yearf <- factor(Data_Geostat$Year)
    Q_ik <- model.matrix(~yearf, contrasts=list(yearf='contr.sum'))
    ## Zero out non-BTS rows so they are unaffected by the lambdas
    Q_ik[which(Data_Geostat$Gear!='Trawl'),] <- 0
  } else {
    ## Constant over time
    Q_ik <- matrix(ifelse(Data_Geostat$Gear=='Trawl', 1, 0), ncol=1)
  }
} else if(model=='ats'){
  ## For this one sum across the two strata to create a single one, akin to
  ## what they'd do without the BTS
  Data_Geostat <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                             Catch_KG=ats$strata2+ats$strata3, depth=ats$depth,
                             depth2=ats$depth2,
                             Gear='Acoustic_3-surface', AreaSwept_km2=1,
                             Vessel='none')
  c_iz <- rep(0, nrow(Data_Geostat))
} else if(model=='bts'){
  Data_Geostat <- DF1
  c_iz <- rep(0, nrow(Data_Geostat))
}
years <- min(Data_Geostat$Year):max(Data_Geostat$Year)
nyr <- length(years)

### Derived objects for spatio-temporal estimation
silent.fn(Spatial_List  <-
  make_spatial_info(grid_size_km=grid_size_km, n_x=n_x,
                    Method=Method, Lon=Data_Geostat[,'Lon'],
                    Lat=Data_Geostat[,'Lat'],
                    fine_scale=finescale,
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
                                   Cov_ep = Data_Geostat[,c('depth', 'depth2')[1]],
                                   Extrapolation_List = Extrapolation_List,
                                   Spatial_List = Spatial_List, FUN = mean,
                                   na.omit = "time-average")))
## Normalize depth and then add depth^2
XX$Cov_xtp <- NULL#(XX$Cov_xtp- mean(XX$Cov_xtp))/sd(XX$Cov_xtp)
## new  <- XX$Cov_xtp[,,1]^2
## XX$Cov_xtp <- abind(XX$Cov_xtp, new, along=3)

## Build data and object for first time
message('Building first TMB object..')
silent.fn(TmbData <- make_data(Version=Version, FieldConfig=FieldConfig,
                   OverdispersionConfig=OverdispersionConfig,
                   RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                   b_i=Data_Geostat[,'Catch_KG'],
                   a_i=Data_Geostat[,'AreaSwept_km2'],
                   ## v_i=as.numeric(Data_Geostat[,'Vessel'])-1,
                   v_i=1:nrow(Data_Geostat),
                   s_i=Data_Geostat[,'knot_i']-1,
                   t_i=Data_Geostat[,'Year'],
                   MeshList=Spatial_List$MeshList,
                   GridList=Spatial_List$GridList,
                   Q_ik=Q_ik,
                   X_gtp=XX$Cov_xtp,
                   spatial_list=Spatial_List,
                   Method=Spatial_List$Method, Options=Options,
                   Aniso=FALSE))
silent.fn(TmbList0 <- make_model(TmbData=TmbData, RunDir=savedir,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir='models', Random="generate"))

## Tweak the Map based on inputs
message("Updating input Map and Params...")
Map <- TmbList0$Map
Params <- TmbList0$Parameters
if(model=='combined'){
  Params$Beta_mean1_c <- c(1.67, -3.1, -3.5)
  Params$Beta_mean2_c <- c(2.8, 3.5,4.2)
}
if(model=='combined'){
  ## Params$L_beta1_z <- c(.2,.3,.5)
  ## Params$L_beta2_z <- c(.6,.3,1)
  if(fixlambda==1) Map$lambda1_k <- factor(NA)
  if(fixlambda==2) Map$lambda2_k <- factor(NA)
  ## both off
  if(fixlambda==12) {
    Map$lambda1_k <- Map$lambda2_k <- factor(NA)
  }
  if(fixlambda==-1) Map$lambda2_k <- factor(NA *Params$lambda2_k)
  Params$logSigmaM[1:3] <- c(1,1,1)
  ## Assume that the two ATS strata have the same observation error
  Map$logSigmaM <- factor( cbind( c(1,2,2), NA, NA) )
} else {
  ##  Params$L_beta1_z <- Params$L_beta2_z <- .4
}
Params$logkappa1 <- Params$logkappa2 <- -5
if(kappaoff==1){
  message("mapping off kappa1")
  Map$logkappa1 <- factor(NA)
} else if(kappaoff==2){
  message("mapping off kappa2")
  Map$logkappa2 <- factor(NA)
} else if(kappaoff==12){
  message("mapping off kappa1 and kappa2")
  Map$logkappa1 <- Map$logkappa2 <- factor(NA)
}

if(space=='ST' & model=='combined'){
  ## Assume rho is the same for strata but only turn on for AR1
  if(length(Params$Beta_rho1_f)!=3) stop('problem with beta_rho1')
  if(length(Params$Beta_rho2_f)!=3) stop('problem with beta_rho2')
  if(beta1temporal & temporal==4) Map$Beta_rho1_f <- factor(c(1,1,1))
  if(beta2temporal & temporal==4) Map$Beta_rho2_f <- factor(c(1,1,1))
}

## Rebuild with the new mapping stuff
TmbList <- make_model(TmbData=TmbData, RunDir=savedir,
                      Version=Version,  RhoConfig=RhoConfig,
                      loc_x=Spatial_List$loc_x, Method=Method,
                      Param=Params, TmbDir='models',
                      Random='generate', Map=Map)
Obj  <-  TmbList[["Obj"]]
if(silent) trash <-  Obj$env$beSilent()

if(model=='combined'){
  message("Reworking L_xx to fix sign switching")
  which.diag <- function(nrow, ncol){
    ## Returns the vector position of the diagonal elements of a nrow x ncol
    ## matrix. This matches the order when VAST converts a vector of L into a matrix
    ## L.  Thus L_vec[which.diag(L_vec)] will be the diagonal elements.
    counter  <- 1
    out <- NULL
    if(ncol=="IID") return(1:3)
    for(r in 1:nrow){
      for(c in 1:ncol){
        ## Only save index of the diagonals
        if(r==c) out <- c(out,counter)
        if(r>=c) counter <- counter+1
      }
    }
    out
  }
  ## Put broad uniform priors on all parameters
  TmbList$Lower[grep('L_omega|L_epsil', names(TmbList$Lower))] <- -10
  TmbList$Upper[grep('L_omega|L_epsil', names(TmbList$Upper))] <- 10

  ## If using multiple factors set teh diagonals to be positive to prevent
  ## label switching
  TmbList$Lower[grep('L_omega1_z', names(TmbList$Lower))[which.diag(3,n_omega1)]] <- 0
  TmbList$Lower[grep('L_omega2_z', names(TmbList$Lower))[which.diag(3,n_omega2)]] <- 0
  TmbList$Lower[grep('L_epsilon1_z', names(TmbList$Lower))[which.diag(3,n_eps1)]] <- 0
  TmbList$Lower[grep('L_epsilon2_z', names(TmbList$Lower))[which.diag(3,n_eps2)]] <- 0
  ## The beta's are just standard deviations in this case so >0
  TmbList$Lower[grep('L_beta1_z', names(TmbList$Lower))] <- 0
  TmbList$Lower[grep('L_beta2_z', names(TmbList$Lower))] <- 0
  ## make sure inits are positive and thus in bound
  par <- Obj$par
  par[grep('L_omega1_z', names(par))[which.diag(3,n_omega1)]]  <-
    abs( par[grep('L_omega1_z', names(par))[which.diag(3,n_omega1)]])
  par[grep('L_omega2_z', names(par))[which.diag(3,n_omega2)]]  <-
    abs( par[grep('L_omega2_z', names(par))[which.diag(3,n_omega2)]])
  par[grep('L_epsilon1_z', names(par))[which.diag(3,n_eps1)]]  <-
    abs( par[grep('L_epsilon1_z', names(par))[which.diag(3,n_eps1)]])
  par[grep('L_epsilon2_z', names(par))[which.diag(3,n_eps2)]]  <-
    abs( par[grep('L_epsilon2_z', names(par))[which.diag(3,n_eps2)]])
  par[grep('L_beta1_z', names(par))] <- abs(par[grep('L_beta1_z', names(par))])
  par[grep('L_beta2_z', names(par))] <- abs(par[grep('L_beta2_z', names(par))])
} else {
  ## These are standard deviations so bound below
  if(!is.na(TmbList$Lower['L_omega1_z'] )) TmbList$Lower['L_omega1_z'] <- 0
  if(!is.na(TmbList$Lower['L_omega2_z'] )) TmbList$Lower['L_omega2_z'] <- 0
  if(!is.na(TmbList$Lower['L_epsilon1_z'] )) TmbList$Lower['L_epsilon1_z'] <- 0
  if(!is.na(TmbList$Lower['L_epsilon2_z'] )) TmbList$Lower['L_epsilon2_z'] <- 0
  if(!is.na(TmbList$Lower['L_beta1_z'] )) TmbList$Lower['L_beta1_z'] <- 0
  if(!is.na(TmbList$Lower['L_beta2_z'] )) TmbList$Lower['L_beta2_z'] <- 0
  ## make sure inits are positive and thus in bound
  par <- Obj$par
  ind <- grep("L_", x=names(par))
  if(length(ind)>0) par[ind] <- abs(par[ind])
}

if(temporal==4){
  TmbList$Upper[grep('rho', names(TmbList$Upper))] <- 1.0
  TmbList$Lower[grep('rho', names(TmbList$Lower))]  <- -1.0
}
## Run it once to optimize the random effects and set that to
## last.par.best which is the init in tmbstan.
Obj$par <- par
Obj$fn(Obj$par)
Obj$env$last.par.best <- Obj$env$last.par

## bundle together some of the inputs that will be needed later for
## plotting and such that aren't included in the standard VAST output
loc <- data.frame(Spatial_List$MeshList$isotropic_mesh$loc[,-3])
names(loc) <- c('E_km', 'N_km')
Inputs <- list(loc=loc, loc_x=data.frame(knot_x=1:n_x, Spatial_List$loc_x))
## MapDetails_list for saving and plotting
mdl <- make_map_info(Region=Region, spatial_list=Spatial_List,
                     Extrapolation_List=Extrapolation_List )
mdl$Legend$x <- mdl$Legend$x-70
mdl$Legend$y <- mdl$Legend$y-45


if(make_plots){
  if(!dir.exists(paste0(savedir, '/data_plots')))
    dir.create(paste0(savedir, '/data_plots'))
  silent.fn(plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List,
                      Data_Geostat=Data_Geostat, PlotDir=paste0(savedir,"/") ))
  ## Some custom maps of the data properties
  if(model=='combined'){
    ## Plot log average catch in grid
    Year_Set <- sort(unique(Data_Geostat$Year))
    Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
    MatDat <- log(tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                         FUN=mean, na.rm=TRUE))
    MatDatSD <- tapply(log(Data_Geostat$Catch_KG), Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                       FUN=sd, na.rm=TRUE)
    ## Some grids have only zero observations
    MatDat[is.infinite(MatDat)]  <-  NA
    MatDatSD[is.infinite(MatDatSD) | is.nan(MatDatSD)]  <-  NA
    ## Use consistent zlim for all three data types
    zlim <- range(MatDat, na.rm=TRUE)
    message('Making data maps by gear type...')
    for(ii in 1:3){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatDat[,ii,Years2Include,drop=TRUE],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/data_plots/map_data_avg_', ii),
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
                 FileName=paste0(savedir, '/data_plots/map_data_pres_', ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=c(0,1),
                 mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
                 textmargin='Presence', zone=mdl$Zone, mar=c(0,0,2,0),
                 oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
    }
    ## Plot standard deviation of data
    for(ii in 1:3){
      PlotMap_Fn(MappingDetails=mdl$MappingDetails,
                 Mat=MatDatSD[,ii,Years2Include,drop=TRUE],
                 PlotDF=mdl$PlotDF,
                 MapSizeRatio=mdl$MapSizeRatio, Xlim=mdl$Xlim, Ylim=mdl$Ylim,
                 FileName=paste0(savedir, '/data_plots/map_data_sd_', ii),
                 Year_Set=Year_Set[Years2Include],
                 Legend=mdl$Legend, zlim=range(MatDatSD, na.rm=TRUE),
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
               FileName=paste0(savedir, '/data_plots/map_data_ratio'),
               Year_Set=Year_Set[Years2Include],
               Legend=mdl$Legend,
               mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include))))),
               textmargin='Ratio log(ATS)/log(BTS)) avg catches', zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE, pch=16)
  }
  ## Calculate raw indices from the data. Took the vAST code and modified to
  ## do it in R. First the totally naive way without space which includes
  ## the added zeroes
  a_g <- as.numeric(TmbData$a_gl)
  IndexNaive <- ddply(subset(Data_Geostat, X != -999), .(Gear,Year), summarize,
                density=sum(a_g)/1000*mean(Catch_KG, na.rm=TRUE))
  gears <- levels(IndexNaive$Gear) ##c('BTS', 'ATS_3-16m', 'ATS>16m')
  D_gcy <- tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                  FUN=mean, na.rm=TRUE)
  Index_cy <- matrix(0, nrow=3, ncol=nyr,
                     dimnames=list(gear=gears,
                                   year=years))
  Index_gcy <- array(NA, dim=c(length(a_g), 3, nyr),
                     dimnames=list(knot=1:length(a_g), gear=gears,
                                   year=years))
  for(y in 1:length(years)){
    ## Expand by area and convert from kg to metric tonnes
    for(cc in 1:3){
      for(g in 1:length(a_g)){
        if(!is.na(D_gcy[g,cc,y])){
          Index_gcy[g,cc,y] <- D_gcy[g,cc,y]*a_g[g]/1000
          Index_cy[cc,y] <- Index_cy[cc,y]+ D_gcy[g,cc,y]*a_g[g]/1000
        }
      }
    }
  }
  index.data <- melt(Index_cy)
  index.data <- index.data[which(index.data$value>0),]
  names(IndexNaive) <- names(index.data)
  index.raw <- rbind(cbind(type='Naive Spatial',index.data),
                     cbind(type='Naive',IndexNaive))
  write.csv(index.raw, file=paste0(savedir, '/index.raw.csv'))
  g <- ggplot(index.raw, aes(year, log(value), group=gear, color=gear)) +
    geom_line() + geom_point() + ylab("Log density") +
    ggtitle(paste0('Raw Data Index w/ n_x=', control$n_x)) + facet_wrap('type') + theme_bw()
  ggsave(paste0(savedir, '/data_plots/raw_data_index.png'), g, width=9, height=5)
  index.data.knot <- dcast(melt(Index_gcy), year+knot~gear)
  png(paste0(savedir, '/data_plots/raw_data_pairs.png'), width=7, height=5,
      units='in', res=500)
  pairs(log(index.data.knot[, c(3:5)]), upper.panel=NULL, cex=.75,
        col=rgb(0,0,0,.5), main=paste('n_x=', control$n_x))
  dev.off()
  ## Also make pairs of presence
  P_gcy <- tapply(Data_Geostat$Catch_KG, Data_Geostat[, c( 'knot_i', 'Gear','Year')],
                  FUN=function(x) mean(x>0, na.rm=TRUE))
  index.presence.knot <- dcast(melt(P_gcy), Year+knot_i~Gear)
  png(paste0(savedir, '/data_plots/raw_data_pairs_presence.png'), width=7, height=5,
      units='in', res=500)
  pairs(index.presence.knot[, c(3:5)], upper.panel=NULL, cex=.75,
        col=rgb(0,0,0,.5), main=paste('n_x=', control$n_x), xlim=c(0,1), ylim=c(0,1))
  dev.off()
}

Record <- list(Version=Version, Method=Method,
               grid_size_km=grid_size_km, n_x=n_x,
               FieldConfig=FieldConfig, RhoConfig=RhoConfig,
               OverdispersionConfig=OverdispersionConfig,
               ObsModel=ObsModel, Region=Region,
               strata.limits=strata.limits,
               control=control, Extrapolation_List=Extrapolation_List,
               Spatial_List=Spatial_List, Data_Geostat=Data_Geostat,
               MapDetails_List=mdl)
save( Record, file=file.path(savedir,"Record.RData"))



### old experimental stuff
## Params$beta2_ft <- Params$beta2_ft+5
## if(model=='combined'){
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
