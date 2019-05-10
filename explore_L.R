
## I need to constrain the L vector elements to avoid sign switching which
## is a problem with MCMC. Not clear how to do this so I'm using simulation
## to explore it. Essentially a prior predictive distribution on the
## covariance matrix.

## What are the implied priors on the covariances and correlations given an
## Lvec?
L_to_covcor <- function(Lvec, nrow=3){
  ## Take L vector and create lower triangular matrix with correlations and
  ## diagonals are variances
  Lvec <- as.numeric(Lvec)
  if(length(Lvec)==5) ncol=2
  if(length(Lvec)==6) ncol=3
  L <- matrix(0, nrow, ncol)
  counter <- 1
  for(r in 1:nrow){
    for(c in 1:ncol){
      if(r>=c){
        L[r,c] <- Lvec[counter]
        counter <- counter+1
      }
    }
  }
  covar <- L %*% t(L)
  covcor <- cov2cor(covar)
  diag(covcor) <- diag(covar)
  covcor ## covcor b/c it's correlation and covariance mixed
}
plot.Leffect <- function(covcor5, covcor6, L5, L6=L5, pch='.'){
  col1 <- rgb(0,0,0,.25)
  col2 <- rgb(1,0,0,.25)
  myplot <- function(i,j, ll, cc) {
    y1 <- if(i==j) sqrt(covcor6[i,j,]) else covcor6[i,j,]
    y2 <- if(i==j) sqrt(covcor5[i,j,]) else covcor5[i,j,]
    plot(x=L6[,ll], y=y1, pch=pch, col=col1, axes=FALSE, ann=FALSE,
         xlim=range(c(L5[,ll], L6[,ll])), ylim=range(c(y1,y2)))
    if(ll!=6){
      points(x=L5[,ll], y=y2, pch=pch, col=col2)
    }
    if(ll==1) axis(2, col=gray(.5))
    if(ll==1) mtext(labs[cc], side=2, line=2)
    if(cc==6) axis(1, col=gray(.5))
    box(col=gray(.5))
  }
  labs <- c("SD(1)", "SD(2)", "SD(3)", "Cor(1,2)", "Cor(1,3)", "Cor(2,3)")
  nL <- ncol(L5)
  par(mfcol=c(6,7), mar=c(0,0,0,0), oma=c(2,4,2,0), mgp=c(1,.5,0))
  ## Add relationships between L and outputs
  for(ll in 1:nL){
    myplot(1,1, ll, 1)
    mtext(paste0("L",ll), line=0)
    myplot(2,2, ll, 2)
    myplot(3,3, ll, 3)
    myplot(1,2, ll, 4)
    myplot(1,3, ll, 5)
    myplot(2,3, ll, 6)
  }
  ## Compare densities of outputs
  myplot2 <- function(i,j, cc){
    if(cc <= 3){
      d5 <- density(covcor5[i,j,], from=0)
      d6 <- density(covcor6[i,j,], from=0)
    } else {
      d5 <- density(covcor5[i,j,], from=-1, to=1)
      d6 <- density(covcor6[i,j,], from=-1, to=1)
    }
    xlim <- range(c(d5$x, d6$x))
    ylim <- c(0, max(c(d5$y, d6$y)))
    plot(0,0, type='n', axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
    box(col=gray(.5))
    lines(d5$x, d5$y, col=col2, lwd=2)
    lines(d6$x, d6$y, col=col1, lwd=2)
  }
  myplot2(1,1,  1)
  mtext('Density')
  myplot2(2,2,  2)
  myplot2(3,3,  3)
  myplot2(1,2,  4)
  myplot2(1,3,  5)
  myplot2(2,3,  6)
}

library(plyr)

ff <- function(file){png(file, width=7, height=5, units='in', res=500)}
nsim <- 15000
set.seed(10)
L <- matrix(runif(nsim*6, min=-10, max=10), ncol=6)
covcor5 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,-6])), perm=c(2,3,1))
covcor6 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,])), perm=c(2,3,1))
ff('plots/Leffect_U10.png')
plot.Leffect(covcor5, covcor6, L)
dev.off()

## Repeat with the bounds currently used
## All Ls are bounded U(-10,10) and some are >0
L[,c(1,3,6)] <- abs(L[,c(1,3,6)])
covcor6 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,])), perm=c(2,3,1))
covcor5 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,-6])), perm=c(2,3,1))
ff('plots/Leffect_diagpos.png')
plot.Leffect(covcor5, covcor6, L)
dev.off()

L[,c(1,3,5,6)] <- abs(L[,c(1,3,5,6)])
covcor6 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,])), perm=c(2,3,1))
covcor5 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,-6])), perm=c(2,3,1))
ff('plots/Leffect_used.png')
plot.Leffect(covcor5, covcor6, L)
dev.off()

## ## Repeat with Normal prior N(0,1)
## L <- matrix(rnorm(nsim*6, 0, sd=1), ncol=6)
## covcor5 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,-6])), perm=c(2,3,1))
## covcor6 <- aperm(laply(1:nsim, function(i) L_to_covcor(L[i,])), perm=c(2,3,1))
## ff('plots/Leffect_N01.png')
## plot.Leffect(covcor5, covcor6, L)
## dev.off()

### Take posterior output and make similar plots
fit <- readRDS('Y:/mcmc_kappaoff_tvlambda_ST/mcmcfit.RDS')
pars <- names(df)[grep('L_omega', names(df))]
png('plots/pairs_temp.png', width=12, height=8, res=800, units='in')
pairs(fit, pars=pars, gap=0)
dev.off()
df <- as.data.frame(fit)
nn <- nrow(df)
set.seed(10)
Lprior <- matrix(runif(nn*6, min=-5, max=5), ncol=6)
Lprior[,c(1,3,6)] <- abs(Lprior[,c(1,3,6)])
Lprior[,6] <- 0

Lpost <- as.matrix(df[, grep('L_omega2_z', x=names(df))])
## Add a fake 6th column so I can use the plotting function.
Lpost <- cbind(Lpost, 0)
covcorprior <- aperm(laply(1:nn, function(i) L_to_covcor(Lprior[i,])), perm=c(2,3,1))
covcorpost <- aperm(laply(1:nn, function(i) L_to_covcor(Lpost[i,])), perm=c(2,3,1))

ff('plots/Leffect_prior_posterior.png')
plot.Leffect(covcorprior, covcorpost, Lprior, Lpost)
dev.off()

## So there's strong bimodality in L and it carries through to Cor(1,2)
## mostly. What is happening here? Let's try making some plots for these
## two regimes.

df1 <- df[df[, grep('L_omega2_z\\[2\\]', names(df))]>0,]
df2 <- df[df[, grep('L_omega2_z\\[2\\]', names(df))]<0,]

source("startup.R")
## Need to manually go to fit_mcmc so the Obj is in the workspace
savedir <- 'L5positive'
dir.create('L5positive')
index <- calculate.index.mcmc(Obj, df1)
index$savedir <- 'L5positive'
plot.index.mcmc(index, 'L5positive')
plot.availability.map.mcmc(index)
plot.density.map.mcmc(index)
plot.covcor.mcmc(index)

## Need to manually go to fit_mcmc so the Obj is in the workspace
savedir <- 'L5negative'
dir.create('L5negative')
index <- calculate.index.mcmc(Obj, df2)
index$savedir <- 'L5negative'
plot.index.mcmc(index, 'L5negative')
plot.availability.map.mcmc(index)
plot.density.map.mcmc(index)
plot.covcor.mcmc(index)

