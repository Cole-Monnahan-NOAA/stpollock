
## Generate a Poisson-link density that matches VAST structure to test for
## consistency
generate.pl.samples <- function(st.list, atrend, nyrs, plot){
  lon <- st.list$lon; lat <- st.list$lat; depth <- st.list$depth; beta0 <- st.list$beta0

  ## long format data.frame to hold the true densities
  D <- data.frame(Lon=lon, Lat=lat, depth=depth)
  D <- D[rep(seq_len(nrow(D)), nyrs), ]
  D$Year <- rep(1:nyrs, each=length(lon))

  ## The Poisson-link predictors
  a1 <- -.15/20; a2 <- .1/20; a3 <- .05/20
  t <- 1:nyrs-1
  beta0 <- 0
  beta1 <- cbind(beta0+t*a1,  beta0+t*a2,  beta0+t*a3)
  ## matplot(beta1)
  beta2 <- matrix(5, nyrs, 3)
  p1 <- p2 <- r1 <- r2 <- den <- matrix(NA, nrow=nrow(D), ncol=3)
  for(i in 1:nrow(D)){
    for(ss in 1:3){
      p1[i,ss] <- beta1[D$Year[i], ss]
      p2[i,ss] <- beta2[ss] + rnorm(1, mean=0, sd=st.list$sd.process)
      r1[i,ss] <- 1-exp(-exp(p1[i,ss]))
      r2[i,ss] <- exp(p1[i,ss]+p2[i,ss])/r1[i,ss]
      den[i,ss] <- exp(p1[i,ss]+p2[i,ss])#r1[i,ss]*r2[i,ss]
    }
  }
  ## matplot(p1)
  ## matplot(p2)
  ## matplot(r1)
  ## matplot(r2)
  ## matplot(den)
  D <- cbind(D, d1=den[,1], d2=den[,2], d3=den[,3])
  D.long <- melt(D, measure.vars=c('d1', 'd2', 'd3'),
                 variable.name='strata', value.name='density')
  if(plot){
    tmp <- ddply(D.long, .(Year, strata), summarize, mean.den=median(density))
    g <- ggplot(tmp, aes(Year, mean.den, color=strata)) + geom_line()
    ggsave(file.path(st.list$plotdir, paste0('pl_annual_density_', st.list$replicate,'.png')), g, width=7, height=5)
  }

  ## now sample from the truth
  i=1
  r1bts <- r1ats1 <- r1ats2 <- dbts <- dats1 <- dats2 <- rep(NA, nrow(D))
  bts <- ats1 <- ats2 <- rep(NA, nrow(D))
  rpl <- function(p, mu, sd){
    ## observed; p is probability of occurence
    x <- rbinom(n=1, size=1, prob=p)
    if(x==1){
      x <- exp(rnorm(1, mean=log(mu), sd=sd)-sd^2/2)
    }
    return(x)
  }
  ## mean(sapply(1:5000, function(i) rpl(.9999, 100, 1.1)))

  for(i in 1:nrow(D)){
    r1bts[i] <- 1-exp(-exp(p1[i,1])-exp(p1[i,2]))
    dbts[i] <- exp(p1[i,1]+p2[i,1])+exp(p1[i,2]+p2[i,2])
    bts[i] <- rpl(r1bts[i], dbts[i]/r1bts[i], sd=st.list$bt.sd)
    r1ats1[i] <- r1[i,2]
    dats1[i] <- exp(p1[i,2]+p2[i,2])
    ats1[i] <- rpl(r1ats1[i], dats1[i]/r1ats1[i], sd=st.list$at.sd)
    r1ats2[i] <- r1[i,3]
    dats2[i] <- exp(p1[i,3]+p2[i,3])
    ats2[i] <- rpl(r1ats2[i], dats2[i]/r1ats2[i], sd=st.list$at.sd)
  }
  D <- cbind(D,  BT=bts, AT1=ats1, AT2=ats2, d1=den[,1], d2=den[,2],
             d3=den[,3])
  pars <- list(beta1_tc=beta1, beta2_tc=beta2, sigma_bt=st.list$bt.sd,
               sigma_at=st.list$at.sd)
  return(list(D=D, pars=pars))
}


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
  ## den2d <- generate.density(st.list=st.list, atrend=atrend,
  ##                           nyrs=st.list$nyrs, plot=plot)
  ## ## Distribution density vertically
  ## den3d <-
  ##   distribute.density(dat=den2d, vtrend=vtrend, st.list=st.list, plot=plot)
  ## ## plot(as.numeric(den3d[1, 6:50]), type='n', ylim=c(0,.5))
  ## ## lapply(sample(1:nrow(den3d), size=10), function(i) lines(as.numeric(den3d[i, 6:50])))

  ## ## Simulate the sampling process for both gear types
  ## data <- sample.data(dat=den3d, st.list=st.list, plot=plot)
  tmp <- generate.pl.samples(st.list=st.list, atrend=atrend, nyrs=st.list$nyrs, plot=plot)
  data <- tmp$D
  truth <- ddply(data, .(Year), summarize,
                 strata1=(sum(d1)),
                 strata2=(sum(d2)),
                 strata3=(sum(d3)))
  ## only makes sense to compare the combiend model I think
  x1 <- data.frame(year=1:st.list$nyrs, par.name='beta1', tmp$pars$beta1_tc)
  x2 <- data.frame(year=1:st.list$nyrs, par.name='beta2',
                   tmp$pars$beta2_tc)
  beta.names <- c('year', 'par.name', 'strata1', 'strata2', 'strata3')
  names(x1) <- names(x2)  <- beta.names
  betas <- melt(rbind(x1,x2), id.vars=c('year', 'par.name'),
                variable.name='strata', value.name='truth')

  ## Fit combinations of models
  k <- 1
  betas.list <- ind.list <- list()
  for(s in spaces){
    for(m in models){
      out <- fit.models(data, replicate, model=m, space=s, plot=plot)
      if(!is.null(out)){
        ## Grab the truth values and the predicted densities by strata to
        ## compare more directly. Have to manually massage the output to
        ## compare to truth
        if(m=='ats') t0 <-log(truth$strata2+truth$strata3)
        if(m=='bts') t0 <-log(truth$strata1+truth$strata2)
        if(m=='combined'){
          b1 <- out$Report$beta1_tc
          b2 <- out$Report$beta2_tc
          x1 <- data.frame(year=1:st.list$nyrs, par.name='beta1', out$Report$beta1_tc)
          x2 <- data.frame(year=1:st.list$nyrs, par.name='beta2', out$Report$beta2_tc)
          names(x1) <- names(x2)  <- beta.names
          betas2 <- melt(rbind(x1,x2), id.vars=c('year', 'par.name'),
                         variable.name='strata', value.name='est')
          stopifnot(identical(betas[,1:3], betas2[,1:3]))
          betas.list[[k]] <- cbind(rep=replicate, model=m, space=s, betas, est=betas2$est)
          tmp <- melt(truth,
                      measure.vars=c('strata1', 'strata2', 'strata3'))
          t0 <- log(tmp$value)
        }
        ind.list[[k]] <-
          cbind(rep=replicate, out$Index, strata2=out$Index.strata$strata,
                strata.est=out$Index.strata$est,
                strata.se=out$Index.strata$se, truth=t0)
      }
      k <- k+1
    }
  }
  ## Return them in long format
  indices <- do.call(rbind, ind.list)
  indices$group <- with(indices, paste(rep, model, strata, sep='-'))
  indices$group2 <- with(indices, paste(rep, model, strata2, sep='-'))
  betas <- do.call(rbind, betas.list)
  betas$group <- with(betas, paste(rep, model, strata, sep='-'))
  return(list(indices=indices, betas=betas))
}



## Basic simulation
set.seed(1)
nyr <- 12
p3 <- seq(.1, .4, len=nyr/2)
p3 <- c(p3, rev(p3))
p2 <- seq(.3, .5, len=nyr/2)
p2 <- c(p2, rev(p2))
p1 <- 1-p3-p2
plot(years, p1, type='n', ylim=c(0,1), xlab=NA, lty=1,
     lwd=2, ylab='Proportion Abundance')
yy <- c(years, rev(years))
polygon(yy, c(rep(0, len=nyr), rev(p1)), col=gray(.2), border=1)
polygon(yy, c(p1,  rev(p1+p2)), col=gray(.5), border=1)
polygon(yy, c(p1+p2,  rev(p1+p2+p3)), col=gray(.8), border=1)
box(col=gray(.5))

atrend <- c(seq(0,1, len=5), seq(1,-1, len=5))
vtrend <- c(4,16,24,16, 20, 25,14, 12,12,12)


## currently depth has no impact but should add that and other covariates later
Nsamples <- 300
st.list <-
  list(lon=runif(Nsamples,-175, -160), lat=runif(Nsamples, 55,62),
       beta0=5, depth = sample(50:51, size=Nsamples, replace=TRUE),
       sd.process=.05, nyrs=10,  bt.sd=.01, at.sd=.01, n_x=10)

## ## Run a single replicate in serial for testing
out <- simulate(replicate=12, st.list=st.list,
                atrend=atrend, vtrend=vtrend, plot=TRUE)
## ggplot(out, aes(year, est, color=strata)) + geom_line() + geom_point() +
##   facet_grid(model~space)

cores <- 6
reps <- cores*2
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
out.parallel <-
  ##sapply(1:reps, function(i)
  snowfall::sfLapply(1:reps, function(i)
    simulate(replicate=i, st.list=st.list, atrend=atrend,
             vtrend=vtrend, plot=TRUE))
sfStop()
indices <- do.call(rbind, lapply(out.parallel, function(i) i$indices))
betas <- do.call(rbind, lapply(out.parallel, function(i) i$betas))

ggplot(indices, aes(year, est, group=group, color=model)) + geom_line()+
  facet_grid(strata~space)
## have to manipulate this one a bit
x <- droplevels(subset(indices, model == 'combined'))
y <- melt(x, id.vars=c('group2', 'year', 'strata2'), measure.vars=c('strata.est', 'truth'))
y$group3 <- paste0(y$group2, '_', y$variable)
ggplot(y, aes(year, value, group=group3, color=variable)) + geom_line() +
  facet_wrap('strata2')

betas2 <- melt(betas, measure.vars=c('truth', 'est'))
betas2$group2 <- paste0(betas2$group, '_', betas2$variable)
ggplot(betas2, aes(year, value, group=group2, color=variable)) + geom_line() +
  facet_grid(par.name~strata, scales='free_y')
