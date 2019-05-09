

bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')

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
## saveRDS(ats.zeroes, file='data/ats.zeroes.RDS')
message("Adding zeroes onto ATS data set")
ats.zeroes <- readRDS('data/ats.zeroes.RDS')
ats <- rbind(ats.zeroes, ats)
## g <- ggplot(ats, aes(lon, lat, col=X==-999)) + geom_point(size=.5) +
##   facet_wrap('year') + theme_bw()
## ggsave('plots/ats.zeroes.png', g, width=10, height=7)

## The ats data is really high resolution so truncating this for now to
## make things faster and fix the mesh which is overly weighted to the ats
## data otherwise
if(filterdata){
  warning("still subsetting the ATS")
  ats <- ats[seq(1, nrow(ats), len=3*nrow(bts)),]
}

## Normalize the covariates. Does it make sense to do that here with depths
## in different strata??
norm <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
## bts$depth <- norm(bts$depth)
bts$depth2 <- bts$depth^2
ats$depth <- ats$surface# norm(ats$surface)
ats$depth2 <- ats$depth^2
ats$mlength <- ats$temp.bottom <- ats$tmp.surface <- NA


## Filter data off the shelf
## hist(bts$depth)
## hist(ats$depth)
## bts <- subset(bts, depth<200)
## ats <- subset(ats, depth<400)

## Temporarily drop some years
## if(filterdata){
##   message("filtering out years <2011")
##   bts <- subset(bts, year <2011)
##   ats <- subset(ats, year <2011)
## }

DF1 <- data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none', depth=bts$depth, depth2=bts$depth2, X=bts$X)
DF2 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, depth2=ats$depth2, X=ats$X)
DF3 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, depth2=ats$depth2, X=ats$X)

## ## Simulate a fake process
## f <- function(mu, sd, p){
##   o <- rbinom(1,1,p)
##   o <- 1
##   if(!o){
##     return(0)
##   } else
##     return(exp(rnorm(1, mu,sd)))
## }
## fake <- data.frame(lat=ats$lat[1:1000], lon=ats$lon[1:1000] , year=rep(1:10, times=100),
##                    strata1=sapply(1:1000, function(i) f(5, .1, .8)),
##                    strata2=sapply(1:1000, function(i) f(6, .1, .8)),
##                    strata3=sapply(1:1000, function(i) f(4, .1, .8)),
##                    depth=0, depth2=0)


## DF1 <- data.frame(Lat=fake$lat, Lon=fake$lon, Year=fake$year,
##                   Catch_KG=fake$strata1+fake$strata2,
##                   Gear='Trawl', AreaSwept_km2=1,
##                   Vessel='none', depth=fake$depth, depth2=fake$depth2)
## DF2 <- data.frame(Lat=fake$lat, Lon=fake$lon, Year=fake$year,
##                   Catch_KG=fake$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
##                   Vessel='none', depth=fake$depth, depth2=fake$depth2)
## DF3 <- data.frame( Lat=fake$lat, Lon=fake$lon, Year=fake$year,
##                    Catch_KG=fake$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
##                    Vessel='none', depth=fake$depth, depth2=fake$depth2)

## DF1$Catch_KG[sample(1:1000, size=80)] <- 0
## DF2$Catch_KG[sample(1:1000, size=80)] <- 0
## DF3$Catch_KG[sample(1:1000, size=80)] <- 0
