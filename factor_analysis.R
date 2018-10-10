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
h1 <- 5; h2 <- 7
ntows <- 355
## Get
at1 <- at2<- rep(NA, ntows)
for (i in 1:ntows){
  at1[i] <- sum(SA[i, 3:16])
  at2[i] <- sum(SA[i, 16:ncol(SA)])
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

