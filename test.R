## Tests to get a TMB version of the model working that matches ADMB
library(TMB)

layers <- as.matrix(read.csv('data/layers.csv', header=FALSE, sep=','))
hauls <- read.csv('data/hauls.csv')

SA <- layers; dimnames(SA) <- NULL
h1 <- 5
h2 <- 7

ntows <- 355
## Prepare layers based on the two h values. This is done internally in
## ADMB but is better outside I think.
sum_SA1 <- sum_SA2 <- rep(NA, ntows)
for (i in 1:ntows){
  sum_SA1[i] <- sum(SA[i, 3:h1])
  sum_SA2[i] <- sum(SA[i, 3:h2])
}

## Prepare the TMB inputs
dat <- list(SA=SA, BSA=hauls$pred_sa, BD=hauls$depth,
            sum_SA1=sum_SA1, sum_SA2=sum_SA2)
pars <- list(log_q=2.5, log_a=8.5,
             d2=-2.1, b_BD=0,
             log_c=3.1,
             logSigma=-.3,
             cb_BD=0)

## Run single iteration to see if NLL matches between TMB and ADMB (with
## same initial values)
compile('models/simplified.cpp')
dyn.load(dynlib('models/simplified'))
obj <- MakeADFun(data=dat, para=pars)
obj$fn()
## ADMB NLL
setwd('models/simplified')
system('simplified -nohess -maxfn 0')
setwd('../..')
## OK they match
## ADMB vs TMB solutions
setwd('models/simplified')
system('simplified')
mle <- R2admb::read_pars('simplified')
setwd('../..')
obj$fn(mle$coefficients[1:7])
mle$loglik
opt <- with(obj, nlminb(mle$coefficients[1:7], fn, gr))
opt$par-mle$coefficients[1:7]
opt <- with(obj, nlminb(pars, fn, gr))
opt$par-mle$coefficients[1:7]
opt$objective
## Can't seem to find the MLE??

## Explore the estimability of this model by starting it from random inits
inits.fn <- function()
  list(log_q=runif(1, -5,5), log_a=runif(1,0,10),
       d2=runif(1, -3,3), b_BD=runif(1,-5,5), log_c=runif(1, 0,6),
       logSigma=runif(1, -2,2), cb_BD=runif(1,-5,5))
solve <- function(){
  tmp <- inits.fn()
  obj$env$beSilent()
  opt <- with(obj, nlminb(start=as.numeric(tmp), objective=fn, gradient=gr, control=list(iter.max=200)))
  ## See if ADMB can do it
  setwd('models/simplified')
  file.remove('simplified.par')
  write.table(x=tmp, file='init.pin', col.names=F, row.names=F)
  system('simplified -ainp init.pin')
  if(file.exists('simplified.par')){
    mle <- R2admb::read_pars('simplified')
    admb <- data.frame(software='admb', nll=-mle$loglik, t(mle$coefficients[1:7]))
  } else {
    admb <- NULL
  }
  setwd('../..')
  obj$fn(mle$coefficients[1:7])
  mle$loglik
  tmb <- data.frame(software='tmb', nll=opt$objective, t(opt$par))
  rbind(tmb, admb)
}

set.seed(3)
fits <- do.call(rbind, lapply(1:3, function(x) solve()))



prof <- tmbprofile(obj, name='log_c', parm.range=c(-20,20))
plot(prof[,1], prof[,2])


library(tmbstan)
options(mc.cores = parallel::detectCores() -1)

fit <- tmbstan(obj, init=function(x) obj$par, iter=1000, chains=3, control=list(max_treedepth=12))
launch_shinystan(fit)
