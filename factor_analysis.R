## This file is to try and get a basic factor analysis working on the
## pollock data where categories are the three depth strata. Based on Jim's
## code from his class.

# Load libraries
library(TMB)
library(TMBhelper)

layers <- as.matrix(read.csv('data/layers.csv', header=FALSE, sep=','))
hauls <- read.csv('data/hauls.csv')
SA <- layers; dimnames(SA) <- NULL
ntows <- 355
## From Stan: first two columns are the ADZ so skip them. 3rd is 0.5-0.75m
## and should assume this is 'h2' or 'h' in the paper. More specifically:
## the first 20 columns are vertical layers from 0 to 5m, every 0.25m, rest
## of the columns are in 1m intervals
layer.widths <- c(seq(.25,5, by=.25), 6+1:(ncol(SA)-20))
atstar <- at1 <- at2<- rep(NA, ntows)
EFH <- 14; iEFH <- which(layer.widths==EFH)
for (i in 1:ntows){
  atstar[i] <- SA[i, 3] # layer just above the ADZ (was sum_SA2)
  at1[i] <- sum(SA[i, 3:iEFH]) # ADZ to EFH (was sum_SA1)
  at2[i] <- sum(SA[i, (iEFH+1):ncol(SA)]) # EFH to surface
}
bt <- hauls$pred_sa
Y <- cbind(log(bt), log(at1), log(at2))
pairs(Y)
plot(log(at1), log(hauls$pred_sa))

## First fit the model with Stan's way (Model D)
dat <- list(BSA=hauls$pred_sa, BD=hauls$depth,
            sum_SA1=at1, sum_SA2=atstar)
pars <- list(log_q=2.5, log_a=8.5, d2=-2.1, b_BD=0, log_c=3.1,
             logSigma=-.3, cb_BD=0)
dyn.unload(dynlib('models/simplified'))
compile('models/simplified.cpp')
dyn.load(dynlib('models/simplified'))
obj1 <- MakeADFun(data=dat, para=pars, DLL='simplified')
obj1$env$beSilent()
opt1 <- Optimize(obj=obj1, getsd=TRUE, newtonsteps=1, control=list(trace=1))
## rep <- sdreport(obj1)
## with(rep, cbind(value, sd))
Report1 <- obj1$report()

### Now fit it with a non-spatial factor analysis using the lognormal distribution
## Prepare the TMB inputs
## Compile if necessary
dyn.unload( dynlib(Version) )
Version <- "models/factor_analysis"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,
dat <- list(Y_sp=Y, n_f=2, X_sj=cbind(rep(1, len=ntows),hauls$depth))
pars <- list(beta_jp=matrix(0,nrow=ncol(dat$X_sj),ncol=ncol(dat$Y_sp)),
              Loadings_vec=rep(1,dat$n_f*ncol(dat$Y_sp)-dat$n_f*(dat$n_f-1)/2),
              "Omega_sf"=matrix(0,nrow=ntows,ncol=dat$n_f),
              logsigma=1)
obj2 <- MakeADFun(data=dat, parameters=pars, random="Omega_sf", hessian=FALSE,
                 inner.control=list(maxit=1000), DLL='factor_analysis')
## table(names(Obj$env$last.par))
obj2$env$beSilent()
# Run model
opt2 <- TMBhelper::Optimize( obj=obj2, getsd=TRUE, newtonsteps=1,
                            control=list(trace=1) )
Report2 <- obj2$report()
## rep <- sdreport(Obj)
## with(rep, cbind(par.fixed, sqrt(diag(cov.fixed))))
## yy <- Report2$logdensity_sp
## par(mfrow=c(1,3))
## plot(yy[,1]+yy[,2], dat$Y_sp[,1]); abline(0,1)
## plot(yy[,2], dat$Y_sp[,2]); abline(0,1)
## plot(yy[,3], dat$Y_sp[,3]); abline(0,1)
## cov.est <- Report2$Loadings_pf%*%t(Report2$Loadings_pf)
## cov2cor(cov.est)


## plot(log(obj$report()$BSA_hat), log(dat$BSA)); abline(0,1)
resids1 <- (Y[,1]-Report2$BT_hat)
resids2 <- (Y[,1]-log(ReportXX$BSA_hat))

par(mfrow=c(3,2))
plot(log(ReportXX$BSA_hat), Y[,1], xlab='Predicted BT',
     ylab='Observed BT', main='Stan model D');abline(0,1)
plot(Report2$BT_hat, Y[,1], xlab='Predicted BT',
     ylab='Observed BT', main='Non-spatial factor analysis');abline(0,1)
plot(log(at1), log(obj$report()$d1), xlab='AT 3-15m',
     ylab='Predicted log ADZ density', )
plot(log(at1), yy[,1], xlab='AT 3-15m',
     ylab='Predicted log ADZ density')
jitter <- rnorm(ntows, 0,.1)
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids2),
     col=ifelse(resids2>0, 'black', 'red'), xlab='long', ylab='lat')
plot(hauls$s_long+jitter, hauls$s_lat, cex=abs(resids1),
     col=ifelse(resids1>0, 'black', 'red'), xlab='long', ylab='lat')
