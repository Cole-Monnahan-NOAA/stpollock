## Tests to get a TMB version of the model working that matches ADMB

library(TMB)

layers <- as.matrix(read.csv('../data/layers.csv', header=FALSE, sep=','))
hauls <- read.csv('../data/hauls.csv')

SA <- layers
h1 <- 5
h2 <- 7

ntows <- 355
## Prepare layers based on the two h values
sum_SA1 <- sum_SA2 <- rep(NA, ntows)
for (i in 1:ntows){
  sum_SA1[i] <- sum(SA[i, 3:h1])
  sum_SA2[i] <- sum(SA[i, 3:h2])
}


dat <- list(SA=SA, BSA=hauls$pred_sa, BD=hauls$depth,
            sum_SA1=sum_SA1, sum_SA2=sum_SA2, logSigma=.123)
pars <- list(log_q=2.5, log_a=8.5,
             d2=0, b_BD=0,
             log_c=3.1,
             logSigma=-.3,
             cb_BD=0)

compile('model.cpp')
dyn.load(dynlib('model'))
obj <- MakeADFun(data=dat, para=pars)
opt <- with(obj, nlminb(pars, fn, gr))
opt$obj
opt$par

prof <- tmbprofile(obj, name='log_c', parm.range=c(-20,20))
plot(prof[,1], prof[,2])


library(tmbstan)
options(mc.cores = parallel::detectCores() -1)

fit <- tmbstan(obj, init=function(x) obj$par, iter=1000, chains=3, control=list(max_treedepth=12))
launch_shinystan(fit)
