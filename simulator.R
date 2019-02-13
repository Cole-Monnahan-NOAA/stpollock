#' Generate spatio-temporal index
#'
#' @return A data.frame containing the spatial locations and densities
#'   associated for each year
generate.density <- function(st.list, atrend, nyrs, plot){
  lon <- st.list$lon; lat <- st.list$lat; depth <- st.list$depth; beta0 <- st.list$beta0
  D <- list()
  for(y in 1:nyrs){
    ## Generate true density, including some real zeroes
    den <- exp( beta0 + atrend[y] + rnorm(n=length(lon), sd=st.list$sd.process))
    den <- den*rbinom(n=length(den), size=1, prob=.9)
    ## for each year create a 2d spatial grid of densities
    D[[y]] <- data.frame(Year=y, Lon=lon, Lat=lat, depth=depth,
                         density=den)
  }
  D <- do.call(rbind, D)
  if(plot){
    g <- ggplot(D, aes(Lon, Lat, size=sqrt(density))) + geom_point() +
      facet_wrap('Year')
    ggsave(file.path(st.list$plotdir, paste0('spatial_density_', st.list$replicate,'.png')), g, width=7, height=5)
  }
  return(D)
}

#' Distribute the spatial density in the vertical dimension
#' @param density The output data.frame from generate.density function.
#' @return The density cbinded with the binned densities
distribute.density <- function(dat, vtrend, st.list, plot, eps=.0001){

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
      (1-p)*dnorm(x=x, mean=vtrend[dat$Year[i]], sd=15)
    ## Truncate really low probabilities to identically 0
    ynorm <- y/sum(y)
    ynorm[ynorm<eps] <- 0
    ynorm <- ynorm/sum(ynorm) # renormalize to be probability
    dvert[i,1:length(y)] <- dat$density[i]*ynorm
  }
  ## check we didn't lose density
  if(max(abs(dat$density-apply(dvert, 1, sum)))>.01)
    warning("Lost density when generating vertical distribution")
  dvert <- as.data.frame(dvert)
  if(plot){
    tows <- rep(seq_along(st.list$lon), times=st.list$nyrs)
    tmp <- reshape2::melt(cbind(year=dat$Year, tow=tows, dvert), id.vars=c('year','tow'))
    tmp$vbin <- as.numeric(tmp$variable)
    tmp <- ddply(tmp, .(year, tow), mutate, rel.density=value/sum(value))
    g <- ggplot(subset(tmp, tow<50), aes(vbin, rel.density, group=tow)) + geom_line() +
      facet_wrap('year')
    ggsave(file.path(st.list$plotdir, paste0('vertical_density_', st.list$replicate,'.png')), g, width=7, height=5)
  }
  names(dvert) <- paste0('d', x.all)
  out <- data.frame(dat, dvert)
  return(out)
}

#' Sample from the 3D density surface for the two gear types.
#'
#' @param density The density output from distribute.density.
#' @return A data frame of locations and sampled gear types
sample.data <- function(dat, st.list, plot){
  bt.sd <- st.list$bt.sd
  at.sd <- st.list$at.sd
  ## Bin down the vertical dimension
  X <- dat[,-(1:5)]
  stopifnot(names(X)[1] == 'd1')
  d1 <- rowSums(X[,1:3])                # ADZ
  d2 <- rowSums(X[,4:16])               # ADZ to EFH
  d3 <- rowSums(X[,-(1:16)])            # EFH to surface
  BT <- exp(rnorm(n=length(d1), mean=log(d1+d2), sd=bt.sd))
  AT1 <- exp(rnorm(n=length(d1), mean=log(d2), sd=at.sd))
  AT2 <- exp(rnorm(n=length(d1), mean=log(d3), sd=at.sd))
  out <- cbind(dat[,1:5],  BT, AT1, AT2, d1, d2, d3)
  if(plot){
    out.long <- melt(out, measure.vars=c('density','BT', 'AT1', 'AT2','d1', 'd2', 'd3'))
    ## sum across observations
    tmp <- ddply(subset(out.long, variable %in% c('density', 'd1', 'd2', 'd3')), .(Year, variable), summarize, abundance=sum(value))
    g <- ggplot(tmp, aes(Year, abundance, group=variable, color=variable)) +
      geom_line()
    ggsave(file.path(st.list$plotdir, paste0('annual_density_', st.list$replicate,'.png')), g, width=7, height=5)
    tmp <- ddply(subset(out.long, variable %in% c('BT', 'AT1', 'AT2')), .(Year, variable), summarize, abundance=sum(value))
    g <- ggplot(tmp, aes(Year, abundance, group=variable, color=variable)) +
      geom_line()
    ggsave(file.path(st.list$plotdir, paste0('annual_observed_', st.list$replicate,'.png')), g, width=7, height=5)
  }
  return(out)
}

#' Prepare simulated inputs for the different model types
#'
#' @param data The output from sample.data
#' @param replicate Replicate number, used for parallel
#' @return A list of fitted objects for each model type
fit.models <- function(data, replicate, model, space, plot){
  ## Setup the VAST model
  ## attach these globally so sourcing works
  DF1 <<- data.frame(Lat=data$Lat, Lon=data$Lon, Year=data$Year,
                    Catch_KG=data$BT, Gear='Trawl', AreaSwept_km2=1,
                    Vessel='none', depth=data$depth, depth2=data$depth^2)
  DF2 <<- data.frame( Lat=data$Lat, Lon=data$Lon, Year=data$Year,
                    Catch_KG=data$AT1, Gear='Acoustic_3-16', AreaSwept_km2=1,
                    Vessel='none', depth=data$depth, depth2=data$depth^2)
  DF3 <<- data.frame( Lat=data$Lat, Lon=data$Lon, Year=data$Year,
                    Catch_KG=data$AT2, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                    Vessel='none', depth=data$depth, depth2=data$depth^2)
  simulated.data <<- TRUE
  model <<- model; space <<- space; n_x <<- st.list$n_x
  if(is.null(n_x)) stop("n_x not defined")
  savedir <<- paste0(getwd(), '/simulations/fit_', replicate, "_",  model, "_", space, '_', n_x)
  source('prepare_inputs.R')
  ## Fit the VAST model
  Opt <- tryCatch(
    Optimize(obj=Obj, lower=TmbList$Lower,
             upper=TmbList$Upper,  savedir=savedir,
             newtonsteps=1, control=list(trace=0)),
             error=function(error_message){
               return(NULL)})
  ## TMBhelper::Check_Identifiable(Obj)
  results <- tryCatch(
    process.results(Opt, Obj, Inputs, model, space, savedir),
    error=function(error_message){
      message(cat('Convergence failed:',replicate, model, space))
      return(NULL)})
  if(!is.null(results) & plot) {
      message("Making plots..")
      tryCatch(plot.vastfit(results), error=function(error_message)
      message(cat('VAST plots failed:',replicate, model, space)))
  }
  return(results)
}


#' Run a single replicate of the simulation
#' @param replicate The replicate number for bookkeeping and setting seed.
#' @param st.list A list of geospatial parameters
#' @param nyrs The number of years to run
#' @param atrend The trend in abundance per year
#' @param vtrend A trend in the mean vertical distribution
## @param bt.sd The variance for the BT sampling process.
## @param at.sd The variance for the AT sampling process.
## @param pl.list Other Poisson-link parameters?
## @param obins The observation bins, assuming 0-3, 3-16, 16-surface
## @param vbins The number of vertical bins used to approximate the
##   vertical density. These are scaled automatically to the water column
##   height (thus bin widths vary with depth).
## @param X A matrix of environmental covariates associated for each
##   spatial point
simulate <- function(replicate, st.list, atrend,
                     vtrend, models=c('ats', 'bts', 'combined'),
                     spaces=c('NS', 'S', 'ST')[1],
                     plot=TRUE){
  ## load libraries again in case run in parallel
  library(VAST); library(TMB); library(TMBhelper); library(ggplot2)
  library(plyr); library(reshape2)

  st.list$plotdir <- file.path(getwd(), 'simulations', 'plots')
  st.list$replicate <- replicate
  ## Generate 2D density
  set.seed(replicate)
  den2d <- generate.density(st.list=st.list, atrend=atrend,
                            nyrs=st.list$nyrs, plot=plot)
  ## Distribution density vertically
  den3d <-
    distribute.density(dat=den2d, vtrend=vtrend, st.list=st.list, plot=plot)
  ## plot(as.numeric(den3d[1, 6:50]), type='n', ylim=c(0,.5))
  ## lapply(sample(1:nrow(den3d), size=10), function(i) lines(as.numeric(den3d[i, 6:50])))

  ## Simulate the sampling process for both gear types
  data <- sample.data(dat=den3d, st.list=st.list, plot=plot)

  ## Fit combinations of models
  k <- 1
  ind.list <- list()
  for(s in spaces){
    for(m in models){
      out <- fit.models(data, replicate, model=m, space=s, plot=plot)
      if(!is.null(out)){
        ind.list[[k]] <- cbind(rep=replicate, out$Index)
      } else {
        ind.list[[k]] <- NULL
      }
      k <- k+1
    }
  }
  ## Return them in long format
  indices <- do.call(rbind, ind.list)
  indices$group <- with(indices, paste(rep, model, strata, sep='-'))
  return(indices)
}

