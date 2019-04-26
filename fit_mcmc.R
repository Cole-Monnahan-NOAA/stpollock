## File to run the fits to the real data
chains <- 6
options(mc.cores = chains)
source('startup.R')
model <- 'combined'


## ST1 with estimated kappa1
control <- list(seed=121, beta2temporal=TRUE, n_x=50, n_eps1=2,
                beta1temporal=TRUE, n_eps2=0, combinedoff=FALSE,
                kappaoff=2, temporal=2, fixlambda=2)
savedir <- paste0(getwd(), '/mcmc_kappa1est_fixlambda1_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=14))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## ST1 w/ kappa1 fixed by timevarying catchability
control <- list(seed=121, beta2temporal=TRUE, n_x=50, n_eps1=2,
                beta1temporal=TRUE, n_eps2=0, combinedoff=FALSE,
                kappaoff=12, temporal=2, fixlambda=-1)
savedir <- paste0(getwd(), '/mcmc_kappaoff_tvlambda_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=14))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## ST1 w/ kappa1 estimated and with timevarying catchability
control <- list(seed=121, beta2temporal=TRUE, n_x=50, n_eps1=2,
                beta1temporal=TRUE, n_eps2=0, combinedoff=FALSE,
                kappaoff=2, temporal=2, fixlambda=-1)
savedir <- paste0(getwd(), '/mcmc_kappa1est_tvlambda_ST')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=600, open_progress=FALSE,
               init='last.par.best',
               control=list(max_treedepth=14))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


## Temp code to plot indices and availability by the three scenarios above
x1 <- readRDS('mcmc_kappa1est_fixlambda1_ST/index.mcmc.RDS')
x2 <- readRDS('mcmc_kappaoff_tvlambda_ST/index.mcmc.RDS')
x3 <- readRDS('mcmc_kappa1est_tvlambda_ST/index.mcmc.RDS')

index.gear <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='By Gear', x$index.gear, scenario=x$scenario)))
index.strata <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='By Strata', x$index.strata, scenario=x$scenario)))
avail <- do.call(rbind, lapply(list(x1,x2,x3), function(x)
  cbind(variable='Availability', x$availability, scenario=x$scenario)))

library(cowplot)

g1 <- ggplot(index.gear, aes(year, y=est, color=gear, group=gear, fill=gear)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('log abundance')
g2 <- ggplot(index.strata, aes(year, y=est, color=stratum, group=stratum, fill=stratum)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('log abundance')
g3 <- ggplot(avail, aes(year, y=est, color=gear, group=gear, fill=gear)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_line(lwd=1.5, alpha=.5)+ theme_bw() +
  facet_wrap('scenario')+ ylab('Availability to gear')
g <- plot_grid(g1, g2,g3, ncol = 1)
save_plot("plots/scenario_grid.png", g, base_width=9, base_height=6,
          base_aspect_ratio = 1.3)

### ## Quick exploration of L_omega1
df <- as.data.frame(fit)
xx <- df[,grep('L_omega1', x=names(df))]
## Each row is a covariance matrix.. see if these are different when L4 is
## positive vs negative
pos <- xx[,4]>0
hist(xx[,4])

L_to_cov <- function(Lvec, ncol=2, nrow=3){
  Lvec <- as.numeric(Lvec)
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
  ## Return vector of just the 6 lower triangular elements
  covar[lower.tri(covar, diag=TRUE)]
}

covs <- data.frame(t(apply(xx, 1, L_to_cov)))
names(covs) <- paste0('cov_vec', 1:6)
png('plots/cov_vec_test.png', width=9, height=7, res=500, units='in')
pairs(covs, col=ifelse(pos, 1,2), gap=0, cex=.5)
dev.off()

## If I arbitrarily set L4 to be positive does it affect the covariances?
xx2 <- xx
xx2[,4] <- abs(xx2[,4])
covs <- data.frame(t(apply(xx2, 1, L_to_cov)))
names(covs) <- paste0('cov_vec', 1:6)
png('plots/cov_vec_test2.png', width=9, height=7, res=500, units='in')
pairs(covs, col=ifelse(pos, 1,2), gap=0, cex=.5)
dev.off()
