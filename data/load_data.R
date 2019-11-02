## This file reads in and prepares the survey data for use in
## VAST.  It is meant to be sourced by other functions.

## Main data sets
bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')

## The "inflated" zeroes (see below too)
if(zeroes.case=='basecase'){
  message("Adding zeroes onto ATS data set")
  ats.zeroes <- readRDS('data/ats.zeroes.basecase.RDS')
} else if(zeroes.case=='sensitivity'){
  ats.zeroes <- readRDS('data/ats.zeroes.sensitivity.RDS')
} else {
  stop("Invalid zeroes.case")
}

ats.zeroes <- ats.zeroes[, names(ats)]
ats <- rbind(ats.zeroes, ats)
ats$depth <- ats$surface# norm(ats$surface)
ats$mlength <- ats$temp.bottom <- ats$tmp.surface <- NA

## Temporarily drop some years
if(filteryears){
  message("filtering out years w/o ATS data")
  ## bts <- subset(bts, year <2011)
  ## ats <- subset(ats, year <2011)
  bts <- subset(bts, year %in% unique(ats$year))
  ## Put these in numeric order to prevent bugs in code later. Otherwise
  ## VAST will make predictions in those years without any data.
  bts$year <- as.numeric(as.factor(bts$year))
  ats$year <- as.numeric(as.factor(ats$year))
}

## Try replicating the first years. This is for troubleshooting
## do not use
if(replicateyears){
  message("replicating early years... only use if you know what you're doing")
  ats2 <- subset(ats, year<=2010) %>% mutate(year=year-4)
  ats3 <- mutate(ats2, year=year-4)
  ats <- rbind(ats, ats2, ats3)
  bts2 <- subset(bts, year<=2010) %>% mutate(year=year-4)
  bts3 <- mutate(bts2, year=year-4)
  bts <- rbind(bts, bts2, bts3)
}

## These get used in prepare_inputs.R script and get passed into
## VAST. These are the three data sets (BT, AT2, AT3).
DF1 <- data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='BT', AreaSwept_km2=1,
                   Vessel='none', depth=bts$depth, X=bts$X)
DF2 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='AT2', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, X=ats$X)
DF3 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='AT3', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, X=ats$X)
message("Done loading data...")

#### --------------------------------------------------
#### Old code to get the AT inflated zeroes. Should not need to
#### redo this.
## ## The hard part is getting the coordinates to plot between. I do this by
## ## plotting the annual data on a map and then use the locator() function to
## ## manually select pairs of points which mimic the track of the AT survey
## ## if it had kept going. Then I place points between each pair. This is
## ## totally manual but I can save the locator output to file for
## ## reproducibility.
## ## To reproduce this, uncomment the below lines and run through each year
## ## using the locator function (see help). Build year at a time and then
## ## combine them together and save to file.
## locs <- list()
## locs[[2007]] <- locate.year(2007)
## locs[[2008]] <- locate.year(2008)
## locs[[2009]] <- locate.year(2009)
## locs[[2010]] <- locate.year(2010)
## locs[[2012]] <- locate.year(2012)
## locs[[2014]] <- locate.year(2014)
## locs[[2016]] <- locate.year(2016)
## locs[[2018]] <- locate.year(2018)
## ats.zeroes <- do.call(rbind, locs)
## Also need to get depths for these fake tows so use the NOAA bath data set
## library(marmap)
## tmpmap <- getNOAA.bathy(lon1 = -168, lon2 = max(ats.zeroes$lon)+1,
##                         lat1 = min(ats.zeroes$lat)-1, lat2 = max(ats.zeroes$lat)+1,
##                         resolution = 1)
## plot.bathy(tmpmap, image=TRUE, land=FALSE, n=50,
##            deepest.isobath=-80, shallowest.isobath=0)
## ## a couple of these look to be on an island and thus are
## ## negative so just drop them
## ats.zeroes$surface <- -1*get.depth(tmpmap, x=ats.zeroes$lon, y=ats.zeroes$lat, locator=FALSE)$depth
## ## old code to check that the predicted depths make sense by
## ## predicting at the real AT data... they are identical
## ## suggesting MACE gets depths from this data set.
## ## ats0 <- subset(ats, lon>= -179)
## ## ats0$depthbath <- -1*get.depth(tmpmap, x=ats0$lon, y=ats0$lat, locator=FALSE)$depth
## ## plot(log(ats0$depth), (ats0$depthbath-ats0$depth))
## ## ggplot(ats0, aes(lon, lat, col=log(depthbath), size=log(depthbath))) + geom_point(alpha=.5) +
## ##   scale_colour_gradient2() + facet_wrap('year') + scale_size(range=c(1,3))
## ## bts0 <- subset(bts, lon>= -179)
## ## bts0$depthbath <- -1*get.depth(tmpmap, x=bts0$lon, y=bts0$lat, locator=FALSE)$depth
## ## plot(log(bts0$depth), (bts0$depthbath-bts0$depth)/bts0$depth)
## ## ggplot(bts0, aes(lon, lat, col=log(depthbath), size=log(depthbath))) + geom_point(alpha=.5) +
## ##   scale_colour_gradient2() + facet_wrap('year') + scale_size(range=c(1,3))
## ## saveRDS(ats.zeroes, file='data/ats.zeroes.RDS')
## ## These are the densities (kg/km2) from Levine & De Robertis
## ## 2019 but previously processed (see manuscript for details).
## inshore <- read.table('data/inshore.catches.txt')[,1]
## ats.zeroes <- readRDS('data/ats.zeroes.full.RDS')
## ats.zeroes$time <- ats.zeroes$date <- ats.zeroes$ground <- NA
## ats.zeroes$strata1 <- NULL
## ## Filter inflated zeroes to be <80m. Then fill in with inshore
## ## observed data from Levine et al 2018 which is from 2017 and
## ## maybe not all pollock so is an upper bound. Thus divide by 4.
## ## plot(lat~lon, data=ats.zeroes, col=ifelse(ats.zeroes$surface>80, 1,2))
## ## ggplot(ats.zeroes, aes(lon, lat, col=surface)) + geom_point()
## ats.zeroes <- subset(ats.zeroes, surface <=80)
## ## Also chop it down to be fewer points
## ats.zeroes <- ats.zeroes[seq(1, nrow(ats.zeroes), by=15),]
## ## I assume that most fish are in the middle stratum, and that
## ## the Levine inshore values are too high so divide those by
## ## 4. Then split the density 4/5 into s2 and 1/5 into s3.
## set.seed(2352152)
## ats.zeroes$strata2 <-
##   rbinom(nrow(ats.zeroes), 1, prob=.4)*sample(inshore*(1/5), size=nrow(ats.zeroes), replace=TRUE)
## ats.zeroes$strata3 <-
##   rbinom(nrow(ats.zeroes), 1, prob=.2)*sample(inshore*(1/20), size=nrow(ats.zeroes), replace=TRUE)
## saveRDS(ats.zeroes, file='data/ats.zeroes.basecase.RDS') # use ats.zeroes.reduced.RDS
## saveRDS(ats.zeroes, file='data/ats.zeroes.sensitivity.RDS') # use ats.zeroes.full.RDS
