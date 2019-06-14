### This files runs the simulation testing component of the analysis.

rm(list=ls())
source("startup.R")
## We condition the OM on the fitted base case model
load('fit_basecase_100/Save.RData')
load('fit_basecase_100/Record.RData')
names(Save)
names(Record)
## Rebuild the Obj with the MLEs
mle <- Save$est$est
par.names <- as.character(Save$est$par)
control <- Record$control
savedir <- paste0(getwd(), '/sim_basecase_100')
source("prepare_inputs.R")
Obj$par <- mle
Obj$fn(mle) -  Save$Opt$objective ## check MLE fit
Obj$gr(mle) ## check that the gradients are 0
pars.all <- Obj$env$last.par ## the full parameter vector

report <- Obj$report(pars.all)
report$Beta_mean1_c
par(mfrow=c(2,2))
matplot(report$beta1_tc)
matplot(report$beta2_tc)
matplot(report$beta1_tc-matrix(report$Beta_mean1_c, ncol=3, nrow=12, byrow=TRUE))
matplot(report$beta2_tc-matrix(report$Beta_mean2_c, ncol=3, nrow=12, byrow=TRUE))


## Put the simulated parameters back into the model
beta1_ft <- (pars.all[grep('beta1_ft', names(pars.all))])
beta1_ft <- matrix(beta1_ft, ncol=3, byrow=TRUE)
as.vector(t(beta1_ft))
par(mfrow=c(1,2))
matplot(beta1_ft)
beta1 <- cbind(c(seq(0,-1.5, len=6), seq(-1.25, 3.5, len=6)),
               c(seq(.5,-1, len=6), seq(-1.25, 2.5, len=6)),
               c(seq(.5,-.5, len=6), seq(-.75,2, len=6)))
matplot(beta1)
pars.all[grep('beta1_ft', names(pars.all))] <- as.vector(t(beta1))

## Simulate the sampling process.
sim <- Obj$report(pars.all)
saveRDS(sim, paste0(savedir, '/simreport.RDS'))
encounters <- rbinom(length(sim$R1_i), size=1, prob=sim$R1_i)
## Carefully index to get the right SigmaObs by gear type
sigma.obs <- sim$SigmaM[,1][as.numeric(Data_Geostat$Gear)]
logcatches <- rnorm(length(sim$R2_i), log(sim$R2_i)-sigma.obs^2/2, sigma.obs)
catches <- exp(logcatches)*encounters


## Rebuild the Obj with the new data
Data_Geostat$Catch_KG <- catches
DF1 <- subset(Data_Geostat, Gear=='Trawl')
DF2 <- subset(Data_Geostat, Gear=='Acoustic_3-16')
DF3 <- subset(Data_Geostat, Gear=='Acoustic_16-surface')
control$simdata <- TRUE
source("prepare_inputs.R")








## Basic simulation
set.seed(1)
nyr <- 12
p3 <- seq(.1, .4, len=nyr/2)
p3 <- c(p3, rev(p3))
p2 <- seq(.3, .5, len=nyr/2)
p2 <- c(p2, rev(p2))
p1 <- 1-p3-p2
plot(years, p1, type='n', ylim=c(0,1), xlab=NA, lty=1,
     lwd=2, ylab='Proportion Abundance')
yy <- c(years, rev(years))
polygon(yy, c(rep(0, len=nyr), rev(p1)), col=gray(.2), border=1)
polygon(yy, c(p1,  rev(p1+p2)), col=gray(.5), border=1)
polygon(yy, c(p1+p2,  rev(p1+p2+p3)), col=gray(.8), border=1)
box(col=gray(.5))

atrend <- c(seq(0,1, len=5), seq(1,-1, len=5))
vtrend <- c(4,16,24,16, 20, 25,14, 12,12,12)


## currently depth has no impact but should add that and other covariates later
Nsamples <- 300
st.list <-
  list(lon=runif(Nsamples,-175, -160), lat=runif(Nsamples, 55,62),
       beta0=5, depth = sample(50:51, size=Nsamples, replace=TRUE),
       sd.process=.05, nyrs=10,  bt.sd=.01, at.sd=.01, n_x=10)

## ## Run a single replicate in serial for testing
out <- simulate(replicate=12, st.list=st.list,
                atrend=atrend, vtrend=vtrend, plot=TRUE)
## ggplot(out, aes(year, est, color=strata)) + geom_line() + geom_point() +
##   facet_grid(model~space)

cores <- 6
reps <- cores*2
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='simulation_progress.txt')
snowfall::sfExportAll()
out.parallel <-
  ##sapply(1:reps, function(i)
  snowfall::sfLapply(1:reps, function(i)
    simulate(replicate=i, st.list=st.list, atrend=atrend,
             vtrend=vtrend, plot=TRUE))
sfStop()
indices <- do.call(rbind, lapply(out.parallel, function(i) i$indices))
betas <- do.call(rbind, lapply(out.parallel, function(i) i$betas))

ggplot(indices, aes(year, est, group=group, color=model)) + geom_line()+
  facet_grid(strata~space)
## have to manipulate this one a bit
x <- droplevels(subset(indices, model == 'combined'))
y <- melt(x, id.vars=c('group2', 'year', 'strata2'), measure.vars=c('strata.est', 'truth'))
y$group3 <- paste0(y$group2, '_', y$variable)
ggplot(y, aes(year, value, group=group3, color=variable)) + geom_line() +
  facet_wrap('strata2')

betas2 <- melt(betas, measure.vars=c('truth', 'est'))
betas2$group2 <- paste0(betas2$group, '_', betas2$variable)
ggplot(betas2, aes(year, value, group=group2, color=variable)) + geom_line() +
  facet_grid(par.name~strata, scales='free_y')
