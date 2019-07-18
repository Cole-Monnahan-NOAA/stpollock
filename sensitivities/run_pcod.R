library(VAST)
library(ggplot2)
library(TMB)

### Try running pcod as a comparable species for prior on logkappa
dat <- FishData::download_catch_rates( survey="EBSBTS", species_set="Gadus macrocephalus" )
dat <- filter(dat, Year>2006)
ggplot(dat, aes(Long, Lat, color=Wt==0)) + geom_point() + facet_wrap('Year')
dir.create('pcod')
setwd('pcod')
settings  <-  make_settings( n_x=100, Region="eastern_bering_sea", purpose="index",
  strata.limits=data.frame(STRATA='All_areas'), bias.correct=FALSE )
settings$RhoConfig[1:4] <- 2
## Run model
fit = fit_model(settings=settings, Lat_i=dat$Lat,
  Lon_i=dat$Long, t_i=dat$Year,
 c_i=rep(0,nrow(dat)),
 b_i=dat$Wt, a_i=rep(1, nrow(dat)))#, "v_i"=dat[,'Vessel'] )
# Plot results
plot( fit )


table <- with(fit$parameter_estimates$SD, data.frame(par=names(par.fixed),
                                                  par.fixed, sd=sqrt(diag(cov.fixed))))
logkappa1 <- as.numeric(table[grep('kappa', table$par),][1,2:3])
logkappa2 <- as.numeric(table[grep('kappa', table$par),][2,2:3])

kappa1 <- exp(1.96*c(1,0,-1)*logkappa1[2]+logkappa1[1])
kappa2 <- exp(1.96*c(1,0,-1)*logkappa2[2]+logkappa2[1])
range1 <- sqrt(8)/kappa1
range2 <- sqrt(8)/kappa2

## > range1
## [1] 334.0540 400.2573 479.5809
## > range2
## [1] 306.5351 351.7551 403.6458
## >

### Use this to calibrate the priors for the H params
rss <- function(V) sqrt(sum(V[1]^2 + V[2]^2))
calculate.H.prior <- function(ln_H, Range1=400, Range2=350){
  H <- matrix(NA, nrow=2, ncol=2)
  H[1,1] = exp(ln_H[1]);
  H[2,1] = ln_H[2];
  H[1,2] = ln_H[2];
  H[2,2] = (1+ln_H[2]*ln_H[2]) / exp(ln_H[1]);
  ## Stol this from plot_anisotropy so could plot all 4 together
  Eigen = eigen(H)
  Pos_Major = Eigen$vectors[, 1] * Eigen$values[1] * Range1
  Pos_Minor = Eigen$vectors[, 2] * Eigen$values[2] * Range1
  Pres_Major = Eigen$vectors[, 1] * Eigen$values[1] * Range2
  Pres_Minor = Eigen$vectors[, 2] * Eigen$values[2] * Range2
  ## Range = 1.1 * c(-1, 1) * max(abs(cbind(Pos_Major, Pos_Minor, Pres_Major, Pres_Minor)))
  ##  print(Range)
  return(list(pos1=Pos_Major, pos2=Pos_Minor, pres1=Pres_Major, pres2=Pres_Minor))
}

draw.prior <- function() rnorm(2, mean= c(0,-1), sd=c(.75,.75))
out <- lapply(1:50, function(i) calculate.H.prior(draw.prior()))
pos1 <- do.call(rbind, lapply(out, function(x) x$pos1))
pos2 <- do.call(rbind, lapply(out, function(x) x$pos2))
pres1 <- do.call(rbind, lapply(out, function(x) x$pres1))
pres2 <- do.call(rbind, lapply(out, function(x) x$pres2))
xlim <- c(-2000, 2000)# range(c(pos1, pres1))
ylim <- c(-2000, 2000)# range(c(pos2, pres2))
pcod.H <- calculate.H.prior(fit$ParHat$ln_H_input)
library(shape)
par(mfrow=c(1,2))
plot(1, type = "n", xlim=xlim , main='Presence', ylim=ylim, xlab = "", ylab = "")
for(i in 1:nrow(pos1)){
plotellipse(rx = rss(pres1[i,]), ry = rss(pres2[i,]),
                   angle = -1*(atan(pres1[i,1]/pres1[i,2])/(2*pi)*360-90),
                   lcol = rgb(0,0,0,.1),
                   lty = c("solid", "dotted")[1])
}
plotellipse(rx = rss(pcod.H$pres1), ry = rss(pcod.H$pres2),
                   angle = -1*(atan(pcod.H$pres1[1]/pcod.H$pres1[2])/(2*pi)*360-90),
                   lcol = 'red',
                   lty = c("solid", "dotted")[1])
plot(1, type = "n", xlim=xlim , ylim=ylim, xlab = "", ylab = "", main='Postive')
for(i in 1:nrow(pos1)){
  plotellipse(rx = rss(pos1[i,]), ry = rss(pos2[i,]),
                     angle = -1 * (atan(pos1[i,1][1]/pos1[i,2])/(2 * pi)
                       * 360 - 90), lcol = rgb(0,0,0,.1),
                     lty = "solid")
}
plotellipse(rx = rss(pcod.H$pos1), ry = rss(pcod.H$pos2),
                   angle = -1*(atan(pcod.H$pos1[1]/pcod.H$pos1[2])/(2*pi)*360-90),
                   lcol = 'red',
                   lty = c("solid", "dotted")[1])
