#' Generate spatio-temporal index
#'
#' @return A data.frame containing the spatial locations and densities
#'   associated for each year
generate.density <- function(st.list, abundance.trend, nyrs, X.space, beta.space){
  lon <- st.list$lon; lat <- st.list$lat; depth <- st.list$depth; beta0 <- st.list$beta0
  D <- list()
  for(y in 1:nyrs){
    ## Generate true density, including some real zeroes
    den <- exp( beta0 + abundance.trend[y] + rnorm(n=length(lon)))
    den <- den*rbinom(n=length(den), size=1, prob=.9)
    ## for each year create a 2d spatial grid of densities
    D[[y]] <- data.frame(Year=y, Lon=lon, Lat=lat, depth=depth,
                         density=den)
  }
  D <- do.call(rbind, D)
  ##  ggplot(D, aes(Lon, Lat, size=sqrt(density))) + geom_point() + facet_wrap('Year')
  return(D)
}

#' Distribute the spatial density in the vertical dimension
#' @param density The output data.frame from generate.density function.
#' @return The density cbinded with the binned densities
distribute.density <- function(dat, vertical.trend, X, vbins, obins,
                               eps=.0001){

  max.depth <- max(dat$depth)
  x.all <- 1:max.depth
  dvert <- matrix(0, nrow=nrow(dat), ncol=max.depth)
  for(i in 1:nrow(dat)){
    ## For each point (row) distribute the density vertically from the
    ## bottom to the surface, which for now is in 1m bins. For now using a
    ## mixture of normals to concentrate fish near bottom
    x <- 1:dat$depth[i] # the vertical bins
    p <- runif(1, .01,.99)
    y <- p*dnorm(x=x, 0, 3) +
      (1-p)*dnorm(x=x, mean=vertical.trend[dat$Year[i]], sd=15)
    ## Truncate really low probabilities to identically 0
    ynorm <- y/sum(y)
    ynorm[ynorm<eps] <- 0
    ynorm <- ynorm/sum(ynorm) # renormalize to be probability
    dvert[i,1:length(y)] <- dat$density[i]*ynorm
  }
  ## check we didn't lose density
  if(max(abs(dat$density-apply(dvert, 1, sum)))>.01)
    warning("Lost density when generating vertical distribution")
  ## quick visual check of distributions
  ##  matplot(t(dvert[1:50,]), type='l')
  dvert <- as.data.frame(dvert)
  names(dvert) <- paste0('d', x.all)
  out <- data.frame(dat, dvert)
  return(out)
}

#' Sample from the 3D density surface for the two gear types.
#'
#' @param density The density output from distribute.density.
#' @return A data frame of locations and sampled gear types
sample.data <- function(dat, bt.sd, at.sd, pl.list, obins){
  ## Bin down the vertical dimension
  X <- dat[,-(1:5)]
  d1 <- rowSums(X[,1:3])                # ADZ
  d2 <- rowSums(X[,4:16])               # ADZ to EFH
  d3 <- rowSums(X[,-(1:16)])            # EFH to surface
  BT <- exp(rnorm(n=length(d1), mean=log(d1+d2), sd=bt.sd))
  AT1 <- exp(rnorm(n=length(d1), mean=log(d2), sd=at.sd))
  AT2 <- exp(rnorm(n=length(d1), mean=log(d3), sd=at.sd))
  out <- cbind(dat[,1:5], BT, AT1, AT2, d1, d2, d3)
  if(FALSE){
    par(mfrow=c(1,3))
    plot(log(d1+d2), log(BT))
    plot(log(d2), log(AT1))
    plot(log(d3), log(AT2))
    out.long <- melt(out, measure.vars=c('d1', 'd2', 'd3'))
    ggplot(out.long, aes(Lon, Lat, size=sqrt(value), col=value==0)) +
      geom_point() + facet_grid(Year~variable)
  }
  return(out)
}

#' Prepare simulated inputs for the different model types
#'
#' @param data The output from sample.data
#' @param replicate Replicate number, used for parallel
#' @return A list of fitted objects for each model type
fit.models <- function(data, replicate, plot=TRUE){

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
  Data_Geostat <-
    reshape2::melt(data, id.vars=c('Lat', 'Lon', 'Year', 'density', 'depth', 'd1', 'd2', 'd3'),
                   value.name='Catch_KG', variable.name='Gear')
  Data_Geostat$Vessel <- factor(1)
  Data_Geostat$AreaSwept_km2 <- 1
  ## Derived objects for spatio-temporal estimation
  DateFile <- paste0(getwd(),'/VAST_output_',replicate, '/')
  Spatial_List  <-  make_spatial_info( grid_size_km=grid_size_km, n_x=n_x,
                                   Method=Method, Lon=Data_Geostat[,'Lon'],
                                   Lat=Data_Geostat[,'Lat'],
                                   Extrapolation_List=Extrapolation_List,
                                   DirPath=DateFile, Save_Results=FALSE )
  Data_Geostat$knot_i=Spatial_List$knot_i
  dir.create(DateFile)
  Record <- list("Version"=Version,"Method"=Method,"grid_size_km"=grid_size_km,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Region"=Region,"strata.limits"=strata.limits)
  save( Record, file=file.path(DateFile,"Record.RData"))
  capture.output( Record, file=paste0(DateFile,"Record.txt"))
  TmbDir <- DateFile
  c_iz = matrix( c(1,2, 2,NA, 3,NA), byrow=TRUE, nrow=3,
                ncol=2)[as.numeric(Data_Geostat[,'Gear']),] - 1
  ## Add threshold
  b_i = Data_Geostat[,'Catch_KG']
  Random = "generate"
  TmbData <- Data_Fn(Version=Version, FieldConfig=FieldConfig,
                     OverdispersionConfig=OverdispersionConfig,
                     RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                     b_i=b_i, a_i=Data_Geostat[,'AreaSwept_km2'],
                     ## v_i=rep(factor(1:110-1), times=3),
                     v_i=Data_Geostat$vessel-1,
                     s_i=Data_Geostat[,'knot_i']-1,
                     t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                     MeshList=Spatial_List$MeshList,
                     GridList=Spatial_List$GridList,
                     Method=Spatial_List$Method, Options=Options,
                     Aniso=FALSE)
  TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                           Version=Version,  RhoConfig=RhoConfig,
                           loc_x=Spatial_List$loc_x, Method=Method,
                           TmbDir='models', Random="generate")
  ## Extract default values
  Map <- TmbList0$Map
  Params <- TmbList0$Parameters
  ## Fix SigmaM for all surveys to be equal
  Map$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
  ##  Map$beta1_ct <- factor(rep(1, 30))
  Map$beta2_ct <- factor(rep(1, 30))
  ## Map$logkappa1 <- factor(NA)
  ## Params$logkappa1 <- 5
  TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                          Version=Version,  RhoConfig=RhoConfig,
                          loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random='generate', Map=Map)

  ## Fit the full VAST model
  obj.full <- TmbList$Obj; obj.full$env$beSilent()
  ## Not sure why passing lower and upper throws an error for this case
  Opt.full <- TMBhelper::Optimize( obj=obj.full, savedir=DateFile, getsd=TRUE,
                   control=list(trace=1))
  rep.full <- obj.full$report()
  fit.full <- list(Opt=Opt.full, Report=rep.full,
                   ## ParHat=obj.full$env$parList(Opt.full$par),
                   TmbData=TmbData)
  if(plot){
    message("Making plots..")
    plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List,
              Data_Geostat=Data_Geostat, PlotDir=DateFile )
    Enc_prob  <-
      plot_encounter_diagnostic(Report=rep.full, Data_Geostat=Data_Geostat,
                                DirName=DateFile)
    Q  <-  plot_quantile_diagnostic( TmbData=TmbData, Report=rep.full, FileName_PP="Posterior_Predictive",
                                    FileName_Phist="Posterior_Predictive-Histogram",
                                    FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile)
    ## Spatial residuals
    MapDetails_List  <- make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
    ## Decide which years to plot
    Year_Set  <-  seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
    Years2Include <-  which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
    plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'],
                   TmbData=TmbData, Report=rep.full, Q=Q, savedir=DateFile,
                   MappingDetails=MapDetails_List[["MappingDetails"]],
                   PlotDF=MapDetails_List[["PlotDF"]],
                   MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
                   Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]],
                   FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include,
                   Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]],
                   Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]],
                   mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)
    Index <-
      plot_biomass_index( DirName=DateFile, TmbData=TmbData,
                         Sdreport=Opt.full[["SD"]], Year_Set=Year_Set,
                         Years2Include=Years2Include,
                         strata_names=strata.limits[,1], use_biascorr=TRUE,
                         category_names=levels(Data_Geostat[,'Gear']) )
    Index$Table[,c("Category","Year","Estimate_metric_tons","SD_mt")]
    Plot_factors( Report=rep.full, ParHat=obj.full$env$parList(), Data=TmbData,
                 SD=Opt$SD, mapdetails_list=MapDetails_List, Year_Set=Year_Set,
                 category_names=levels(Data_Geostat[,'Gear']), plotdir=DateFile )
  }

  ## Now two independent ST indices with VAST


  ## And the combined model from Stan's paper
  ##  dyn.unload(dynlib(paste0(DateFile, 'VAST_v4_0_0')))
  return(list(vast.full=fit.full))
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
#' @param bt.sd The variance for the BT sampling process.
#' @param at.sd The variance for the AT sampling process.
#' @param pl.list Other Poisson-link parameters?
#' @param obins The observation bins, assuming 0-3, 3-16, 16-surface
simulate <- function(replicate, st.list, nyrs, abundance.trend,
                     vertical.trend, vbins, X, bt.sd, at.sd, pl.list,
                     obins){
  ## load libraries again in case run in parallel
  library(VAST); library(TMB); library(TMBhelper)
  ## Check inputs TODO

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

  bt.sd <- .2
  at.sd <- .2
  data <- sample.data(den3d, bt.sd, at.sd, pl.list, obins)

  ## Reorganize data for models
  out <- fit.models(data, replicate, plot=TRUE)




  ## Fit the independent VAST models, one for BT and one for AT

  ## and the combined model

  return(out)

}

