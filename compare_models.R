## This file is to try and get a basic factor analysis working on the
## pollock data where categories are the three depth strata. Based on Jim's
## code from his class.

## Load and prepare data for TMB
n_f <- 3 # number of factors
source("load_data.R")

## First fit the model with Stan's way (Model D)
dat <- list(BSA=Y[,1], BD=hauls$depth, sum_SA1=at1, sum_SA2=atstar)
pars <- list(log_q=2.5, log_a=8.5, d2=-2.1, b_BD=0, log_c=3.1,
             logSigma=-.3, cb_BD=0)
compile('models/simplified.cpp')
dyn.load(dynlib('models/simplified'))
obj1 <- MakeADFun(data=dat, para=pars, DLL='simplified')
obj1$env$beSilent()
opt1 <- Optimize(obj=obj1, getsd=TRUE, newtonsteps=1, control=list(trace=0))
Report1 <- obj1$report()

### Now fit it with a non-spatial factor analysis using the lognormal
### distribution
Version <- "models/factor_analysis"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
dat <- list(Y_sp=Y, n_f=n_f, X_sj=cbind(rep(1, len=ntows),hauls$depth))
pars <- list(beta_jp=matrix(0,nrow=ncol(dat$X_sj),ncol=ncol(dat$Y_sp)),
              Loadings_vec=lvec,
              "Omega_sf"=matrix(0,nrow=ntows,ncol=dat$n_f),
              logsigma=1)
obj2 <- MakeADFun(data=dat, parameters=pars, random="Omega_sf",
                  DLL='factor_analysis')
obj2$env$beSilent()
opt2 <- Optimize( obj=obj2, getsd=TRUE, newtonsteps=1,
                            control=list(trace=0) )
Report2 <- obj2$report()

### Now fit it with a non-spatial factor analysis using the Poisson-link
### likelihood
Version <- "models/factor_analysis_pois"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
dat <- list(Y_sp=(Y), n_f=n_f, X_sj=cbind(rep(1, len=ntows),hauls$depth))
pars <- list(beta_jp=matrix(.1,nrow=ncol(dat$X_sj),ncol=ncol(dat$Y_sp)),
              Loadings_vec=lvec,
              "Omega_sf"=matrix(0,nrow=ntows,ncol=dat$n_f),
              logsigma=1, logweight=-6)
obj3 <- MakeADFun(data=dat, parameters=pars, random="Omega_sf",
                  DLL='factor_analysis_pois',
                  map=list(logweight=factor(NA)))
obj3$env$beSilent()
opt3 <- Optimize( obj=obj3, getsd=TRUE, control=list(trace=0))
Report3 <- obj3$report()

### Now fit it with a spatial factor analysis using the Poisson-link
### likelihood
Version <- "models/spatial_factor_analysis_pois"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
dat <- list(Y_sp=(Y), n_f=n_f, n_x=mesh$n, x_s=mesh$idx$loc-1,
            X_sj=cbind(rep(1, len=ntows),hauls$depth),
            M0=spde$param.inla$M0, M1=spde$param.inla$M1,
            M2=spde$param.inla$M2)
pars <- list(beta_jp=matrix(.1, nrow=ncol(dat$X_sj),ncol=ncol(dat$Y_sp)),
              Loadings_vec=lvec,
              "Omega_xf"=matrix(0,nrow=dat$n_x,ncol=dat$n_f),
              logsigma=1, logweight=-6, log_kappa=.1)
obj4 <- MakeADFun(data=dat, parameters=pars, random="Omega_xf",
                  DLL='spatial_factor_analysis_pois',
                  map=list(logweight=factor(NA)))
obj4$env$beSilent()
opt4 <- Optimize(obj=obj4, getsd=TRUE, control=list(trace=0))
Report4 <- obj4$report()

## Try turning off the FA part of the SFA model
obj5 <- MakeADFun(data=dat, parameters=pars, random="Omega_xf",
                  DLL='spatial_factor_analysis_pois',
                  map=list(logweight=factor(NA), Loadings_vec=fa.off))
obj5$env$beSilent()
opt5 <- Optimize(obj=obj5, getsd=TRUE, control=list(trace=0))
Report5 <- obj5$report()

## Try turning off the spatial part by setting log_kappa really big but
## turning FA part back on
pars$log_kappa <- 5
obj6 <- MakeADFun(data=dat, parameters=pars, random="Omega_xf",
                  DLL='spatial_factor_analysis_pois',
                  map=list(logweight=factor(NA), log_kappa=factor(NA)))
obj6$env$beSilent()
opt6 <- Optimize(obj=obj6, getsd=TRUE, control=list(trace=0))
opt6 <- opt6$opt # not converging so returned list is different
Report6 <- obj6$report()


## Model comparisons
reps <- list(Report1, Report2, Report3, Report4, Report5, Report6)
tmp <- do.call(rbind, lapply(4:6, function(i)
  cbind(model=i, lon=mesh$loc[,1], lat=mesh$loc[,2], reps[[i]]$Omega_xf)))
sr <- melt(data.frame(tmp), id.vars=c('model', 'lon', 'lat'),
  varnames=c('strata1', 'strata2', 'strata3'))
ggplot(sr, aes(lon, lat, size=abs(value), col=value>0)) +
  geom_point(alpha=.5) +
  facet_grid(model~variable)

fits <- list(opt1, opt2, opt3, opt4, opt5, opt6)
aics <- sapply(fits, function(x) round(x$AIC,2))
npars <- sapply(fits, function(x) x$number_of_coefficients[2])
nlls <- sapply(fits, function(x) round(x$objective,2))
results <- data.frame(rbind(npars, nlls, aics))
mnames <- c('Stan', 'FA:lognormal', 'FA:Poisson', 'SFA:Poisson',
                    'Spatial only', 'FA only')
names(results) <- mnames
results

## These are the predicted densities in the ADZ
y1 <- log(Report1$d1)
y2 <- Report2$logdensity_sp
y3 <- Report3$logdensity
y4 <- Report4$logdensity
## The residuals (log scale)
resids1 <- (Y[,1]-Report1$BSA_hat)
resids2 <- (Y[,1]-Report2$BT_hat)
resids3 <- (Y[,1]-Report3$BT_hat)
resids4 <- (Y[,1]-Report4$BT_hat)


png('plots/method_comparison.png', width=6.5, height=6, res=500, units='in')
par(mfrow=c(3,4), mgp=c(1,.2,0), mar=c(2,3,2,.5), tck=-.02)
xlim <- range(c(Report1$BSA_hat, Report2$BT_hat, Report3$BT_hat))
plot(Report1$BSA_hat, Y[,1], xlab='Predicted BT', xlim=xlim,
     ylab='Observed BT', main='Stan model D');abline(0,1)
plot(Report2$BT_hat, Y[,1], xlab='Predicted BT', xlim=xlim,
     ylab='Observed BT', main='Non-SFA: lognormal');abline(0,1)
plot(Report3$BT_hat, Y[,1], xlab='Predicted BT', xlim=xlim,
     ylab='Observed BT', main='Non-SFA: Poisson-link');abline(0,1)
plot(Report4$BT_hat, Y[,1], xlab='Predicted BT', xlim=xlim,
     ylab='Observed BT', main='SFA: Poisson-link');abline(0,1)
ylim <- range(c(y1, y2[,1], y3[,1], y4[,1]))
plot(log(at1), y1, xlab='Observed log AT 3-15m',
     ylab='Predicted log ADZ density', ylim=ylim)
abline(0,1)
plot(log(at1), y2[,1], xlab='Observed log AT 3-15m',
     ylab='Predicted log ADZ density', ylim=ylim)
abline(0,1)
plot(log(at1), y3[,1], xlab='Observed log AT 3-15m',
     ylab='Predicted log ADZ density', ylim=ylim)
abline(0,1)
plot(log(at1), y4[,1], xlab='Observed log AT 3-15m',
     ylab='Predicted log ADZ density', ylim=ylim)
abline(0,1)
set.seed(20)
jitter <- rnorm(ntows, 0,.1)
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids1),
     col=ifelse(resids1>0, 'black', 'red'), xlab='long', ylab='lat')
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids2),
     col=ifelse(resids2>0, 'black', 'red'), xlab='long', ylab='lat')
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids3),
     col=ifelse(resids3>0, 'black', 'red'), xlab='long', ylab='lat')
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids4),
     col=ifelse(resids4>0, 'black', 'red'), xlab='long', ylab='lat')
dev.off()

png('plots/method_comparison_spatial_resids.png', width=5, height=6, res=500, units='in')
par(mfrow=c(2,2), mgp=c(1,.2,0), mar=c(.5,.5,.5,.5), tck=-.02)
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids1), axes=FALSE,
     col=ifelse(resids1>0, 'black', 'red'), xlab='long', ylab='lat')
mtext('Stan D', line=-1.5); box()
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids2), axes=FALSE,
     col=ifelse(resids2>0, 'black', 'red'), xlab='long', ylab='lat')
mtext('FA: lognormal', line=-1.5); box()
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids3), axes=FALSE,
     col=ifelse(resids3>0, 'black', 'red'), xlab='long', ylab='lat')
mtext('FA: Poisson-link', line=-1.5); box()
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids4), axes=FALSE,
     col=ifelse(resids4>0, 'black', 'red'), xlab='long', ylab='lat')
mtext('SFA: Poisson-link', line=-1.5); box()
dev.off()

## make some plots to explore the spatial patterns
par(mfrow=c(3,3), mgp=c(1,.2,0), mar=c(2,3,2,.5), tck=-.02)
x <- hauls$s_long; y <- hauls$s_lat
Z <- apply(Y, 2, function(x) sqrt(exp(x))/20)
for(i in 1:3) plot(x,y, cex=Z[,i])
## The non-spatial factor analysis
Z <- apply(Report3$logdensity, 2, function(x) sqrt(exp(x))/20)
for(i in 1:3) plot(x,y, cex=Z[,i])
## SFA
Z <- apply(Report4$logdensity, 2, function(x) sqrt(exp(x))/20)
for(i in 1:3) plot(x,y, cex=Z[,i])

## where are the differences from adding space?
par(mfrow=c(1,3))
for(i in 1:3){
  plot(Report3$logdensity[,i], Report4$logdensity[,i])
  abline(0,1)
}



### Now try fitting the same exact data in VAST
TmbList = Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir='models', Random=Random)
# Extract default values
Map = TmbList$Map
Params = TmbList$Parameters
# Fix SigmaM for all surveys to be equal
Map$logSigmaM = factor( cbind( c(1,1,1), NA, NA) )

# Re-build object
TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile,
                       "Version"=Version,  "RhoConfig"=RhoConfig,
                       "loc_x"=Spatial_List$loc_x, "Method"=Method,
                       "TmbDir"=TmbDir, "Random"=Random, Map=Map)
Obj = TmbList[["Obj"]]
Obj$env$beSilent()
## Estimate fixed effects and predict random effects
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]],
                          upper=TmbList[["Upper"]], getsd=TRUE,
                          newtonsteps=1, savedir=DateFile,
                          bias.correct=FALSE ,
                          control=list(trace=1))
ReportVast = Obj$report()
