library(VAST)

### Try running pcod as a comparable species for prior on logkappa
dat <- FishData::download_catch_rates( survey="EBSBTS", species_set="Gadus macrocephalus" )
dat <- filter(dat, Year>2006)
ggplot(dat, aes(Long, Lat, color=Wt==0)) + geom_point() + facet_wrap('Year')

## Make settings (turning off bias.correct to save time for example)
dir.create('pcod')
setwd('pcod')
settings = make_settings( n_x=100, Region="eastern_bering_sea", purpose="index",
  strata.limits=example$strata.limits, bias.correct=FALSE )
settings$RhoConfig[1:4] <- 2
## Run model
fit = fit_model(settings=settings, Lat_i=dat$Lat,
  Lon_i=dat$Lon, t_i=dat$Year,
 c_i=rep(0,nrow(dat)),
 b_i=dat$Wt, a_i=rep(1, nrow(dat)))#, "v_i"=dat[,'Vessel'] )

# Plot results
plot( fit )

logkappa1 <- c(-4.50852914, 0.17220052)
logkappa2 <- c(-4.76112353, 0.14143555)

kappa1 <- exp(1.96*c(1,0,-1)*logkappa1[2]+logkappa1[1])
kappa2 <- exp(1.96*c(1,0,-1)*logkappa2[2]+logkappa2[1])
range1 <- sqrt(8)/kappa1
range2 <- sqrt(8)/kappa2
