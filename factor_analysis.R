## This file is to try and get a basic factor analysis working on the
## pollock data where categories are the three depth strata. Based on Jim's
## code from his class.

# Load libraries
library(TMB)
## Compile if necessary
dyn.unload( dynlib(Version) )
Version = "models/factor_analysis"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

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
Y <- cbind(log(hauls$pred_sa), log(at1), log(at2))
pairs(Y)
## Prepare the TMB inputs
dat <- list(Y_sp=Y, n_f=2, X_sj=matrix(1, nrow=ntows))
Params = list(beta_jp=matrix(0,nrow=ncol(dat$X_sj),ncol=ncol(dat$Y_sp)),
              Loadings_vec=rep(1,dat$n_f*ncol(dat$Y_sp)-dat$n_f*(dat$n_f-1)/2),
              "Omega_sf"=matrix(0,nrow=ntows,ncol=dat$n_f),
              logsigma=1)
# Declare random
Random = c("Omega_sf")

# Initialization
Obj <- MakeADFun(data=dat, parameters=Params, random=Random, hessian=FALSE, inner.control=list(maxit=1000) )
table(names(Obj$env$last.par))
Obj$env$beSilent()
# Run model
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1, control=list(trace=1) )

# Summarize
Report = Obj$report()
rep <- sdreport(Obj)
with(rep, cbind(par.fixed, sqrt(diag(cov.fixed))))

yy <- Report$ln_yexp_sp
par(mfrow=c(1,3))
plot(yy[,1]+yy[,2], dat$Y_sp[,1]); abline(0,1)
plot(yy[,2], dat$Y_sp[,2]); abline(0,1)
plot(yy[,3], dat$Y_sp[,3]); abline(0,1)

cov.est <- Report$Loadings_pf%*%t(Report$Loadings_pf)
cov2cor(cov.est)


#### Fit the same data with Stan's model
dat <- list(BSA=hauls$pred_sa, BD=hauls$depth,
            sum_SA1=at1, sum_SA2=atstar)
pars <- list(log_q=2.5, log_a=8.5,
             d2=-2.1, b_BD=0,
             log_c=3.1,
             logSigma=-.3,
             cb_BD=0)

## Run single iteration to see if NLL matches between TMB and ADMB (with
## same initial values)
dyn.unload(dynlib('models/simplified'))
compile('models/simplified.cpp')
dyn.load(dynlib('models/simplified'))
obj <- MakeADFun(data=dat, para=pars, DLL='simplified')
obj$env$beSilent()
Opt = TMBhelper::Optimize(obj=obj, getsd=TRUE, newtonsteps=1, control=list(trace=1))
rep <- sdreport(obj)
with(rep, cbind(value, sd))
plot(log(obj$report()$d1), log(at1))
plot(yy[,1], log(at1))
