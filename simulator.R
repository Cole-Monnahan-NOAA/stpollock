#' Generate spatio-temporal index
#'
#' @return A data.frame containing the spatial locations and densities
#'   associated for each year
generate.density <- function(st.list, abundance.trend, nyrs, X.space, beta.space){
  lon <- st.list$lon; lat <- st.list$lat; depth <- st.list$depth
  D <- list()
  for(y in 1:nyrs){
    D[[y]] <- data.frame(Year=y, Lon=lon, Lat=lat, depth=depth,
          density=exp(abundance.trend[y] + rnorm(n=length(lon))))
  }
  D <- do.call(rbind, D)
  return(D)
}

#' Distribute the spatial density in the vertical dimension
#' @param density The output data.frame from generate.density function.
#' @return The density cbinded with the binned densities
distribute.density <- function(density, vertical.trend, X, vbins, obins,
                               eps=.0001){
  dvert <- matrix(0, nrow=nrow(density), ncol=max(density$depth))
  for(i in 1:nrow(density)){
    x <- 1:density$depth[i]
    y <- dnorm(x, 0, 1) +
      dnorm(x=x, mean=vertical.trend[density$Year[i]], sd=5)
    ## Truncate really low probabilities to identically 0
    ynorm <- y/sum(y)
    ynorm[ynorm<eps] <- 0
    ynorm <- ynorm/sum(ynorm) # renormalize to be probability
    dvert[i,1:length(y)] <- density$density[i]*ynorm
  }
  ## max(abs(density$density-apply(dvert, 1, sum)))
  out <- data.frame(density, dvert)
  return(out)
}

#' Sample from the 3D density surface for the two gear types.
#'
#' @param density The density output from distribute.density.
#' @return A data frame of locations and sampled gear types
sample.data <- function(density, bt.cv, at.cv, pl.list, obins){
  ## Bin down the vertical dimension
  X <- density[,-(1:5)]
  d1 <- rowSums(X[,1:3])                # ADZ
  d2 <- rowSums(X[,4:16])               # ADZ to EFH
  d3 <- rowSums(X[,-(1:16)])            # EFH to surface
  BT <- exp(rnorm(n=length(d1), mean=log(d1+d2), sd=.2))
  AT1 <- exp(rnorm(n=length(d1), mean=log(d2), sd=.2))
  AT2 <- exp(rnorm(n=length(d1), mean=log(d3), sd=.2))
  out <- cbind(density[,1:5], BT, AT1, AT2)
  return(out)
}

#' Prepare simulated inputs for the different model types
#'
#' @param data The output from sample.data
#' @return A list of lists for the data inputs for all of the model types
prepare.inputs <- function(data){

  ## Setup the VAST inputs
  Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
  grid_size_km <- 50
  n_x <- 100
  ## Stratification for results
  strata.limits <- data.frame(STRATA='All_areas')
  ## Derived objects
  Region <- "Eastern_Bering_Sea"
  ## Save settings
  Extrapolation_List  <-
    make_extrapolation_info( Region=Region, strata.limits=strata.limits )

  ## Model settings
  FieldConfig <- c(Omega1=3, Epsilon1=0, Omega2=0, Epsilon2=0)
  RhoConfig <- c(Beta1=0, Beta2=0, Epsilon1=0, Epsilon2=0)
  OverdispersionConfig <- c(Delta1=0, Delta2=0)
  ObsModel <- c(1,1)
  Options <-  c(SD_site_density=0, SD_site_logdensity=0, Calculate_Range=1,
                Calculate_evenness=0, Calculate_effective_area=1, Calculate_Cov_SE=0,
                Calculate_Synchrony=0, Calculate_Coherence=0)
  Data_Geostat <- melt(data, id.vars=c('Lat', 'Lon', 'Year', 'density', 'depth'),
                       value.name='Catch_KG',
                       variable.name='Gear')
  Data_Geostat$Vessel <- factor(1)
  Data_Geostat$AreaSwept_km2 <- 1
  ## Derived objects for spatio-temporal estimation
  Spatial_List  <-  make_spatial_info( grid_size_km=grid_size_km, n_x=n_x,
                                   Method=Method, Lon=Data_Geostat[,'Lon'],
                                   Lat=Data_Geostat[,'Lat'],
                                   Extrapolation_List=Extrapolation_List,
                                   DirPath=DateFile, Save_Results=FALSE )
  Data_Geostat$knot_i=Spatial_List$knot_i
  DateFile <- paste0(getwd(),'/VAST_output/')
  dir.create(DateFile)
  Record <- list("Version"="VAST_v4_0_0","Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
  save( Record, file=file.path(DateFile,"Record.RData"))
  capture.output( Record, file=paste0(DateFile,"Record.txt"))
  TmbDir <- DateFile

  c_iz = matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
                ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
                                        # Add threshold
  b_i = Data_Geostat[,'Catch_KG']
  Random = "generate"
  TmbData <- Data_Fn(Version="VAST_v4_0_0", FieldConfig=FieldConfig,
                      OverdispersionConfig=OverdispersionConfig,
                      RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                      b_i=b_i, a_i=Data_Geostat[,'AreaSwept_km2'],
                      v_i=rep(factor(1:110-1), times=3),
                      s_i=Data_Geostat[,'knot_i']-1,
                      t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                      MeshList=Spatial_List$MeshList,
                      GridList=Spatial_List$GridList,
                      Method=Spatial_List$Method, Options=Options,
                      Aniso=FALSE)
  TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                           Version="VAST_v4_0_0",  RhoConfig=RhoConfig,
                           loc_x=Spatial_List$loc_x, Method=Method,
                           TmbDir='models', Random="generate")
  ## Extract default values
  Map <- TmbList0$Map
  Params <- TmbList0$Parameters
  ## Fix SigmaM for all surveys to be equal
  Map$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
  TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                          Version="VAST_v4_0_0",  RhoConfig=RhoConfig,
                          loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random='generate', Map=Map)


  ## Now two independent ST indices with VAST


  ## And the combined model from Stan's paper


  return(list(vast.full=TmbList))
}


#' Run a single replicate of the simulation
#' @param replicate The replicate number for bookkeeping and setting seed.
#' @param st.list A list of geospatial parameters
#' @param nyrs The number of years to run
#' @param abundance.trend The trend in abundance per year
#' @param vertical.trend A trend in the mean vertical distribution
#' @param vbins The number of vertical bins used to approximate the
#'   vertical density. These are scaled automatically to the water column
#'   height (thus bin widths vary with depth).
#' @param X A matrix of environmental covariates associated for each
#'   spatial point
#' @param bt.cv The variance for the BT sampling process.
#' @param at.cv The variance for the AT sampling process.
#' @param pl.list Other Poisson-link parameters?
#' @param obins The observation bins, assuming 0-3, 3-16, 16-surface
simulate <- function(replicate, st.list, nyrs, abundance.trend,
                     vertical.trend, vbins, X, bt.cv, at.cv, pl.list,
                     obins){
  ## Check inputs

  ## Generate 2D density
  set.seed(replicate)
  den2d <- generate.density(st.list=st.list, abundance.trend=atrend,
                            nyrs=nyrs, X.space=NULL, beta.space=NULL)

  ## Distribution density vertically
  vertical.trend <- c(4,16,24,16, 18, 18,14, 12,12,12)
  den3d <- distribute.density(den2d, vertical.trend, X.vert, beta.vert,
                              obins)
  ## plot(as.numeric(den3d[1, 6:50]), type='n', ylim=c(0,.5))
  ## lapply(sample(1:nrow(den3d), size=10), function(i) lines(as.numeric(den3d[i, 6:50])))

  ## Simulate the sampling process for both gear types

  bt.cv <- .2
  at.cv <- .2
  data <- sample.data(den3d, bt.cv, at.cv, pl.list, obins)

  ## Reorganize data for models
  inputs <- prepare.inputs(data)

  ## Fit the full VAST model
  ## x <- inputs$vast.full
  ## build.full <- Build_TMB_Fn(TmbData=x$TmbData, RunDir=x$DateFile,
  ##                         Version=x$Version,  RhoConfig=x$RhoConfig,
  ##                         loc_x=x$Spatial_List$loc_x, Method=x$Method,
  ##                         TmbDir=x$TmbDir, Random=x$Random)
  obj.full <- inputs$vast.full$Obj; obj.full$env$beSilent()
  ## Not sure why passing lower and upper throws an error for this case
  Opt.full <- Optimize( obj=obj.full, savedir=DateFile, getsd=TRUE,
                       ##                   lower=build.full$lower, upper=build.full$upper,
                   control=list(trace=0))
  rep.full <- obj.full$report()
  fit.full <- list(index=apply(rep.full$Index_cyl, 2, sum))
  ## Fit the independent VAST models, one for BT and one for AT

  ## and the combined model

  return(list(fit.full))

}


set.seed(1)
atrend <- sort(rnorm(10))
out <- simulate(replicate=1,
                st.list=list(lon=runif(100), lat=runif(100),
                             depth = sample(50:100, size=100, replace=TRUE)),
                nyrs=10, abundance.trend=atrend)
plot(atrend,out[[1]]$index)
