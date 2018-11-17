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

### --------------------------------------------------
### Now try fitting the same exact data in VAST
Version <- "VAST_v4_0_0"
TmbList0 <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir='models', Random=Random)
## Extract default values
Map <- TmbList0$Map
Params <- TmbList0$Parameters
## Fix SigmaM for all surveys to be equal
Map$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random, Map=Map)
Obj1 <- TmbList[["Obj"]]; Obj1$env$beSilent()
Opt1 <- Optimize( obj=Obj1, lower=TmbList[["Lower"]],
                 upper=TmbList[["Upper"]], savedir=DateFile,
                 control=list(trace=1))
ReportVast1 <- Obj1$report()

## Refit but try turning off spatial impact by setting logkappa big
Params$logkappa1 <- 5
Map$logkappa1 <- factor(NA)
TmbList <- Build_TMB_Fn(TmbData=TmbData, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random, Map=Map)
Obj2 <- TmbList[["Obj"]]; Obj2$env$beSilent()
Opt2 <- Optimize( obj=Obj2, lower=TmbList[["Lower"]],
                 upper=TmbList[["Upper"]], savedir=DateFile,
                 control=list(trace=1))
ReportVast2 <- Obj2$report()


## Refit without space at all and putting loadings matrix on overdispersion
FieldConfig <- c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
TmbData3 <- Data_Fn(Version="VAST_v4_0_0", FieldConfig=FieldConfig,
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
## The function breaks below so manually construct the Par and Map inputs
## Setup the new loadings
TmbList <- Build_TMB_Fn(TmbData=TmbData3, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random, build_model=FALSE)
x <- TmbList$Map
y <- TmbList$Parameters
for(i in names(x)){
  if(length(x[[i]]) != length(y[[i]]))
    print(i)
}
## this is a bug in VAST?
x$L2_z <- factor(NA)
x$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
## Rebuild with altered map
TmbList3 <- Build_TMB_Fn(TmbData=TmbData3, RunDir=DateFile,
                       Map=x, Param=y,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random)
Obj3 <- TmbList3[["Obj"]]; Obj3$env$beSilent()
## Not sure why passing lower and upper throws an error for this case
Opt3 <- Optimize( obj=Obj3, savedir=DateFile, getsd=FALSE,
                 control=list(trace=1))
ReportVast3 <- Obj3$report()

## now a version with just space
FieldConfig <- c("Omega1"=3, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig <- c("Delta1"=0, "Delta2"=0)
TmbData4 <- Data_Fn(Version="VAST_v4_0_0", FieldConfig=FieldConfig,
                  OverdispersionConfig=OverdispersionConfig,
                  RhoConfig=RhoConfig, ObsModel=ObsModel, c_iz=c_iz,
                  b_i=b_i, a_i=Data_Geostat[,'AreaSwept_km2'],
                  v_i=as.numeric(Data_Geostat[,'Vessel'])-1,
                  s_i=Data_Geostat[,'knot_i']-1,
                  t_i=Data_Geostat[,'Year'], a_xl=Spatial_List$a_xl,
                  MeshList=Spatial_List$MeshList,
                  GridList=Spatial_List$GridList,
                  Method=Spatial_List$Method, Options=Options,
                  Aniso=FALSE)
## The function breaks below so manually construct the Par and Map inputs
## Setup the new loadings
TmbList0 <- Build_TMB_Fn(TmbData=TmbData4, RunDir=DateFile,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random, build_model=FALSE)
TmbList0$Map$logSigmaM <- factor( cbind( c(1,1,1), NA, NA) )
## Rebuild with altered map
TmbList4 <- Build_TMB_Fn(TmbData=TmbData3, RunDir=DateFile,
                       Map=TmbList0$Map, Param=TmbList0$Parameters,
                       Version=Version,  RhoConfig=RhoConfig,
                       loc_x=Spatial_List$loc_x, Method=Method,
                       TmbDir=TmbDir, Random=Random)
Obj3 <- TmbList3[["Obj"]]; Obj3$env$beSilent()
## Not sure why passing lower and upper throws an error for this case
Opt3 <- Optimize( obj=Obj3, savedir=DateFile, getsd=FALSE,
                 control=list(trace=1))
ReportVast3 <- Obj3$report()


## Look at GRMFs for the different combinations. Have to do some crazy
## stuff to get them to compare to VAST.
xx <- Spatial_List$MeshList$isotropic_mesh$loc
x1 <- data.frame(model='VAST: full', E_km=xx[,1], N_km=xx[,2], ReportVast1$Omegainput1_sf)
x2 <- data.frame(model='VAST: high kappa', E_km=xx[,1], N_km=xx[,2],
                 ReportVast2$Omegainput1_sf)
x3 <- data.frame(model='VAST: no space', E_km=hauls$s_long, N_km=hauls$s_lat,
                 ReportVast3$eta1_vf)
UTMlist <- Convert_LL_to_UTM_Fn( Lon=x3$E_km, Lat=x3$N_km,
                                zone=Extrapolation_List$zone,
                                flip_around_dateline=Extrapolation_List$flip_around_dateline )
x3$E_km <- UTMlist$X
x3$N_km <- UTMlist$Y
xx <- rbind(x1,x2,x3)
reps <- list(Report1, Report2, Report3, Report4, Report5, Report6)
tmp <- data.frame(do.call(rbind, lapply(4:6, function(i)
  cbind(model=i, E_km=mesh$loc[,1], N_km=mesh$loc[,2], reps[[i]]$Omega_xf))))
UTMlist <- Convert_LL_to_UTM_Fn( Lon=tmp[,2], Lat=tmp[,3],
                                zone=Extrapolation_List$zone,
                                flip_around_dateline=Extrapolation_List$flip_around_dateline )
tmp$E_km <- UTMlist$X
tmp$N_km <- UTMlist$Y
mnames <- c('Stan', 'FA:lognormal', 'FA:Poisson', 'SFA:Poisson',
                    'Spatial only', 'FA only')
tmp$model <- factor(mnames[tmp$model])
names(xx) <- names(tmp) <- c("model", 'E_km', 'N_km', 'grmf1', 'grmf2', 'grmf3')
tmp2 <- rbind(xx,tmp)
sr <- melt(data.frame(tmp2), id.vars=c('model', 'E_km', 'N_km'))
ggplot(sr, aes(E_km, N_km, size=abs(value), col=value<0)) +
  geom_point(alpha=.5) +
  facet_grid(model~variable)


fits <- list(opt1, opt2, opt3, opt4, opt5, opt6)
aics <- sapply(fits, function(x) round(x$AIC,2))
npars <- sapply(fits, function(x) x$number_of_coefficients[2])
nlls <- sapply(fits, function(x) round(x$objective,2))
results <- data.frame(rbind(npars, nlls, aics))
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



