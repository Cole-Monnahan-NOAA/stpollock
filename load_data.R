

bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')
set.seed(1)
ats$surface <- rnorm(nrow(ats))
bts$depth <- rnorm(nrow(bts))

## The ats data is really high resolution so truncating this for now to
## make things faster and fix the mesh which is overly weighted to the ats
## data otherwise
ats <- ats[seq(1, nrow(ats), len=nrow(bts)),]

## Normalize the covariates
norm <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
bts$depth <- norm(bts$depth)
bts$depth2 <- bts$depth^2
## bts$mlength <- norm(bts$mean_length)
## bts$temp.bottom <- norm(bts$temp.bottom)
## bts$temp.surface <- norm(bts$temp.surface)
ats$depth <- norm(ats$surface)
ats$depth2 <- ats$depth^2
ats$mlength <- ats$temp.bottom <- ats$tmp.surface <- NA

DF1 <- data.frame( Lat=bts$lat, Lon=bts$lon, Year=bts$year,
                   Catch_KG=bts$density, Gear='Trawl', AreaSwept_km2=1,
                   Vessel='none', depth=bts$depth, depth2=bts$depth2)
DF2 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata2, Gear='Acoustic_3-16', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, depth2=ats$depth2)
DF3 <- data.frame( Lat=ats$lat, Lon=ats$lon, Year=ats$year,
                   Catch_KG=ats$strata3, Gear='Acoustic_16-surface', AreaSwept_km2=1,
                   Vessel='none', depth=ats$depth, depth2=ats$depth2)

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
