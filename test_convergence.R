## script to dig deeper into convergence issues by looking deeper at betas ñ

source("startup.R")
## Test combined spatial model
n_x <- 50
model <- 'combined'; space <- 'ST'
control <- list(seed=112, n_eps2=0, beta2temporal=FALSE,
                finescale=TRUE)
savedir <- paste0(getwd(), '/fit_', model, "_", space,  "_", n_x)
##set.seed(112) ## seed 111 works for ST; 112 crashes out
options(warn=0)
source("prepare_inputs.R")
## options(warn=2) # stop on a warning
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=TRUE, loopnum=10,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=0, control=list(iter.max=300, trace=1))
options(warn=0)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)

## Check where problem occurred
all <- Obj$env$last.par
fixed <- all[-Obj$env$random]
random <- all[Obj$env$random]
Obj$fn(fixed)
Obj$gr(fixed)
jnll.crash <- Obj$report(all)$jnll
hes <- Obj$env$spHess(par=all, random=TRUE)
evs <- eigen(hes)
sort(evs$values)[1:4]
png('testing/eigenvalues.png', res=500, units='in', width=7, height=5)
plot(log10(evs$values+1e-4))
dev.off()

## Try rebuilding it with LA turned off and fixed effects at the MLE
Map$L_beta1_z <- factor(c(NA, NA, NA))
Map$L_epsilon1_z <- factor(rep(NA,5))
Map$lambda1_k <- factor(Params$lambda1_k*NA)
Params$lambda1_k <- fixed['lambda1_k']
Map$L_omega1_z <- factor(Params$L_omega1_z*NA)
Params$L_omega1_z <- fixed[grep('L_omega1_z', x=names(fixed))]
Map$L_epsilon1_z <- factor(Params$L_epsilon1_z*NA)
Params$L_epsilon1_z <- fixed[grep('L_epsilon1_z', x=names(fixed))]
Map$L_beta1_z <- factor(Params$L_beta1_z*NA)
Params$L_beta1_z <- fixed[grep('L_beta1_z', x=names(fixed))]
Map$logkappa1 <- factor(Params$logkappa1*NA)
Params$logkappa1 <- fixed[grep('logkappa1', x=names(fixed))]
Map$Beta_mean1_c <- factor(Params$Beta_mean1_c*NA)
Params$Beta_mean1_c <- fixed[grep('Beta_mean1_c', x=names(fixed))]
Map$Beta_rho1_f <- factor(Params$Beta_rho1_f*NA)
Params$Beta_rho1_f <- fixed[grep('Beta_rho1_f', x=names(fixed))]
Map$lambda2_k <- factor(Params$lambda2_k*NA)
Params$lambda2_k <- fixed[grep('lambda2_k', x=names(fixed))]
Map$L_omega2_z <- factor(Params$L_omega2_z*NA)
Params$L_omega2_z <- fixed[grep('L_omega2_z', x=names(fixed))]
Map$logkappa2 <- factor(Params$logkappa2*NA)
Params$logkappa2 <- fixed[grep('logkappa2', x=names(fixed))]
Params$beta1_ft <- matrix(random[grep('beta1_ft', x=names(random))], nrow=3)
## This one is different
Map$beta2_ft <- factor(as.numeric(TmbList0$Map$beta2_ft) * NA)
Params$beta2_ft <- matrix(fixed[grep('beta2_ft', x=names(fixed))][TmbList0$Map$beta2_ft],nrow=3)
Map$logSigmaM <- factor(Params$logSigmaM*NA)
Params$logSigmaM <- TmbList0$Parameters$logSigmaM
Params$logSigmaM[1,1] <- fixed[grep('logSigmaM', x=names(fixed))][1]
Params$logSigmaM[2:3,1] <- fixed[grep('logSigmaM', x=names(fixed))][2]
Map$Epsilon_rho1_f <- factor(c(NA, NA))
Params$Epsilon_rho1_f <- rep(fixed[grep('Epsilon_rho1_f', x=names(fixed))],2)
Map$Beta_rho1_f <- factor(c(NA, NA,NA))
Params$Beta_rho1_f <- rep(fixed[grep('Beta_rho1_f', x=names(fixed))],3)
Params$Epsiloninput1_sft <-
  array(all[grep('Epsiloninput1_sft', names(all))], dim=c(66,2,12))
Params$Omegainput1_sf <-  matrix(all[grep('Omegainput1_sf', names(all))],
                                 nrow=66, ncol=3)
Params$Omegainput2_sf <-  matrix(all[grep('Omegainput2_sf', names(all))],
                                 nrow=66, ncol=3)
TmbList2 <- make_model(TmbData=TmbData, RunDir=savedir,
                      Version=Version,  RhoConfig=RhoConfig,
                      loc_x=Spatial_List$loc_x, Method=Method,
                      Param=Params, TmbDir='models',
                      Random=NULL, Map=Map)
Obj2  <-  TmbList2[["Obj"]]
Obj2$env$beSilent()
## plot(Obj2$par-random)                   #make sure worked
Obj2$fn()-jnll.crash
grs <- as.numeric(Obj2$gr())
range(grs)

## ## Try optimizing to see if bad Hessian still
## Opt2 <- Optimize(obj=Obj2, lower=TmbList2$Lower, loopnum=3, getsd=FALSE,
##                 upper=TmbList2$Upper,  savedir=savedir,
##                 newtonsteps=0, control=list(trace=20))
bad <- Check_Identifiable(Obj2)$BadParams
bad.pars <- which(bad$Param_check=='Bad')
bad[bad.pars,]
bad.pars.names <- paste(rownames(bad), bad$Param[bad.pars], sep="_")
hes <- Obj2$he(Obj2$par)
evs <- eigen(hes)
sort(evs$values)[1:4]
png('testing/eigenvalues2.png', res=500, units='in', width=7, height=5)
plot(log10(evs$values+1e-4))
dev.off()

##What do they look like spatially?
allpars <- Obj2$env$last.par
omega2 <- array(allpars[grep('Omegainput2', x=names(allpars))],
              dim=c(66,3))
badomega2 <- array('Bad'!=bad$Param_check[grep('Omegainput2', x=names(allpars))],
              dim=c(66,3))
dimnames(badomega2) <- list(knot=1:66, factor=c('factor1', 'factor2', 'factor3'))
temp <- cbind(model='combined', space='ST', knot=1:66, Inputs$loc)
baddf <- merge(melt(badomega2, value.name='Identifiable'), temp,  by='knot')
g <- ggplot(baddf, aes(E_km, N_km, col=Identifiable)) +
  geom_point(size=2) + facet_grid(factor~.)
ggsave(file.path('testing', 'bad_omega2_by_knot.png'), g, width=6, height=4, units='in')

eps1 <- array(allpars[grep('Epsiloninput1', x=names(allpars))],
              dim=c(66,2,12))
badeps1 <- array('Bad'!=bad$Param_check[grep('Epsiloninput1', x=names(allpars))],
              dim=c(66,2,12))
dimnames(badeps1) <- list(knot=1:66, factor=c('factor1', 'factor2'),
                       year=1:12)
temp <- cbind(model='combined', space='ST', knot=1:66, Inputs$loc)
baddf <- merge(melt(badeps1, value.name='Identifiable'), temp,  by='knot')
g <- ggplot(subset(baddf, year>3 & year <10),
            aes(E_km, N_km, col=Identifiable)) +
  geom_point(size=2) +
  facet_grid(year~factor)
ggsave(file.path('testing', 'bad_epsilon1_by_knot.png'), g, width=7, height=7, units='in')

## How about the surface of the bad parameters using MCMC?
library(tmbstan)
library(shinystan)
chains <- 7
options(mc.cores = chains)
inits <- function() Opt2$par
fit2 <- tmbstan(Obj2,  chains=chains,
               iter=800, open_progress=FALSE,
               init=inits, control=list(max_treedepth=10))
saveRDS(object = fit2, file='fit2.RDS')
fit2 <- readRDS('fit2.RDS')
launch_shinystan(fit2)
png('testing/pairs_badpars.png', res=500, units='in', width=12, height=10)
pars <- names(fit2)[bad.pars]
pairs(fit2, pars=pars[1:7+11], gap=0)
dev.off()



## Try using MCMC to see if it helps identify the issue
options(warn=0)
library(tmbstan)
library(shinystan)
## The only tricky part is getting the L's to be identifiable b/c of label
## switching. I think if I set the diagonals to be positive everything else
## will fall into place
lwr <- rep(-Inf, length(all)); upr <- -1*lwr
names(lwr) <- names(upr) <- names(all)
lwr[which(names(all) %in% names(TmbList$Lower))] <- TmbList$Lower
upr[which(names(all) %in% names(TmbList$Upper))] <- TmbList$Upper
## Now the variances
## names(all) <- paste0(1:length(all),"_", names(all))
## names(all)[grep('L_', x=names(all))]
lwr[c(38,40,43)] <- 0 # digonal for L_omega1_z n_f=3
lwr[c(44,46)] <- 0 # digonal for L_epsilon_z which is n_f=2
lwr[c(1844,1846,1849)] <- 0 # diagonal for L_omega2_z n_f=3
lwr[c(49,50,51)] <- 0 # sd for the temporal rw on beta1s
## make sure inits are positive for the random effects, the difference in
## index is b/c the order in all is not the same
all[c(38,40,43,44,46,1844,1846,1849,49,50,51)] <-
  abs(all[c(38,40,43,44,46,1844,1846,1849,49,50,51)])
inits <- function() all
chains <- 7
options(mc.cores = chains)
fit <- tmbstan(Obj, lower=lwr, upper=upr, chains=chains,
               iter=1200, open_progress=FALSE,
               init=inits, control=list(max_treedepth=12))
saveRDS(object = fit, file='fit.RDS')
fit <- readRDS('fit.RDS')
launch_shinystan(fit)

mon <- monitor(fit, print=FALSE)
stanpars <- dimnames(mon)[[1]]
slow <- sort(mon[,'n_eff'])[1:10]
barplot(slow)
png('testing/mcmc_pairs_slow.png', width=12, height=7, units='in', res=500)
slow <- names(sort(mon[,'n_eff'])[1:8])
pairs(fit, pars=slow, gap=0, upper.panel=NULL)
dev.off()
png('testing/mcmc_pairs.png', width=7, height=5, units='in', res=500)
pairs(fit, gap=0, pars=c('logkappa1','logkappa2', 'beta2_ft[1]', 'beta2_ft[2]',
                  'beta2_ft[3]', 'Epsilon_rho1_f', 'Beta_rho1_f', 'lp__'))
dev.off()
for(ii in 1:3){
  png(paste0('testing/mcmc_pairs_beta1_ft_',ii,'.png'), width=15, height=10, units='in', res=500)
  pars <- stanpars[grep('beta1_ft', stanpars)][seq(1,34, by=3)+ii-1]
  pairs(fit, gap=0, pars=pars)
  dev.off()
}

## Try to find a standard model that doesn't work to show Jim.
run.iteration <- function(seed){
  set.seed(seed)
  n_x <<- 50
  model <<- 'combined'; space <<- 'ST'
  combinedoff <<- FALSE
  savedir <<- paste0(getwd(), '/testing/test_std_', seed)
  source('startup.R')
  source("prepare_inputs.R")
  options(warn=2) # stop immediately on NaN warning to save time
  err <- tryCatch(Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=FALSE,
                                  loopnum=3,
                                  upper=TmbList$Upper,  savedir=savedir,
                                  newtonsteps=0, control=list(iter.max=120, trace=1)),
                  error=function(e) NULL)
  if(is.null(err)){
    all <- Obj$env$last.par
    fixed <- all[-Obj$env$random]
    ## Crashed out so save par values to test where
    out <- list(crashed=TRUE, all=all, fixed=fixed, init=Obj$par)
  } else {
    out <- list(crashed=FALSE, max_gradient=Opt$max_gradient,
                objective=Opt$objective, par=Opt$par, init=Obj$par)
  }
  dyn.unload(dynlib(paste0(savedir,'/VAST_v8_0_0')))
  unlink(savedir, recursive=TRUE)
  return(out)
}

library(snowfall)
cores <- 20
chains <- cores*10
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='convergence_progress.txt')
sfExport('run.iteration')
out.parallel <- sfLapply(1:chains, function(i) run.iteration(i))
saveRDS('testing/out.parallel.RDS', object=out.parallel)
sfStop()

out.parallel <- readRDS('testing/out.parallel.RDS')
## First look at converged cases
conv.ind <- which(sapply(out.parallel, function(x) !x$crashed))

png('testing/inits_conv.png', width=7, height=5, units='in', res=500)
par(mfrow=c(1,2))
mgs <- sapply(out.parallel[conv.ind], function(x) log10(x$max_gradient))
mnlls <- sapply(out.parallel[conv.ind], function(x) x$objective)
plot(mgs, ylab='log10 max gradient')
plot(mnlls, ylab='marginal NLL')
dev.off()

parsall <- data.frame(t(sapply(out.parallel, function(x) x$init)))
crashed <- sapply(out.parallel, function(x) x$crashed)
mean(crashed)
## resort to improve plot
ind <- order(crashed, decreasing=TRUE)
parsall <- parsall[ind,]
crashed <- crashed[ind]
## Some aren't different so drop them
parsall2 <- parsall[,-which(apply(parsall, 2, sd)==0)]
png('testing/inits_crashed.png', res=500, units='in', width=7, height=3.5)
par(mar=c(4,4,1,1))
plot(0,0, xlim=c(1, ncol(parsall2)+1), ylim=range(parsall2), type='n',
     xlab='Parameter number', ylab='Initial value')
for(i in 1:nrow(parsall2)){
  points(x=1:ncol(parsall2)+crashed[i]*.25, y=parsall2[i,], col=ifelse(crashed[i], 'black', 'red'))
}
dev.off()


## Took the console trace output and processed it into Excel to read back
## in. Did this for a failed and successful run to compare. This takes like
## half an hour to run (maybe more?)
pars0 <- read.table('testing/trace_crash.csv', sep=',')
grads0 <- jnll0 <- mnll0 <- list()
for(i in 1:nrow(pars0)){
  p  <- as.numeric(pars0[i,])
  mnll0[[i]] <-  as.numeric(Obj$fn(p))
  grads0[[i]] <- Obj$gr(p)
  jnll0[[i]] <- Obj$report(Obj$env$last.par)$jnll
  print(i)
}
mnll0 <- do.call(c, mnll0)
grads0 <- as.data.frame(do.call(rbind, grads0))
jnll0 <- do.call(c, jnll0)
names(grads0) <- names(pars0) <- paste0(1:length(Obj$par),"_", names(Obj$par))
maxgrads0 <- apply(grads0, 1, function(x) max(abs(x)))
grads0$iter <- pars0$iter <- 1:nrow(grads0)
pars0.long <- melt(pars0, id.vars='iter', value.name='value')
grads0.long <- melt(grads0, id.vars='iter', value.name='gradient')
pars1 <- read.table('testing/trace_nocrash.csv', sep=',')
grads1 <- jnll1 <- mnll1 <- list()
for(i in 1:nrow(pars1)){
  p  <- as.numeric(pars1[i,])
  mnll1[[i]] <-  as.numeric(Obj$fn(p))
  grads1[[i]] <- Obj$gr(p)
  jnll1[[i]] <- Obj$report(Obj$env$last.par)$jnll
  print(i)
}
mnll1 <- do.call(c, mnll1)
grads1 <- as.data.frame(do.call(rbind, grads1))
jnll1 <- do.call(c, jnll1)
names(grads1) <- names(pars1) <- paste0(1:length(Obj$par),"_", names(Obj$par))
maxgrads1 <- apply(grads1, 1, function(x) max(abs(x)))
grads1$iter <- pars1$iter <- 1:nrow(grads1)
pars1.long <- melt(pars1, id.vars='iter', value.name='value')
grads1.long <- melt(grads1, id.vars='iter', value.name='gradient')
## Now combine them together and make plots
pars.long <- rbind(cbind(converged=FALSE, pars0.long),
                   cbind(converged=TRUE, pars1.long))
g <- ggplot(pars.long, aes(iter, value, color=converged)) + geom_line(lwd=2,alpha=.5) +
                   facet_wrap('variable', scales='free_y') + theme_bw()
ggsave('testing/pars_crash_nocrash.png', g, width=12, height=7)
grads.long <- rbind(cbind(converged=FALSE, grads0.long),
                   cbind(converged=TRUE, grads1.long))
g <- ggplot(grads.long, aes(iter, gradient, color=converged)) + geom_line(lwd=2,alpha=.5) +
                   facet_wrap('variable', scales='free_y') + theme_bw() +
  geom_abline(slope=0, intercept=0)
ggsave('testing/grads_crash_nocrash.png', g, width=12, height=7)
likes <- as.data.frame(rbind(data.frame(iter=1:length(jnll0), converged=FALSE, joint=jnll0, marginal=mnll0),
               data.frame(iter=1:length(jnll1), converged=TRUE, joint=jnll1, marginal=mnll1)))
likes.long <- melt(likes, id.vars=c('iter', 'converged'),
               value.name='negloglike', variable.name='type')
g <- ggplot(likes.long, aes(iter, negloglike, color=converged))+ geom_line(lwd=2) +
  facet_wrap('type', scales='free_y', ncol=1)
ggsave('testing/nlls_crash_nocrash.png', g, width=7, height=5)



## Look at data by knot more closely
tmp <- ddply(Data_Geostat, .(knot_i, Gear, Year), summarize,
             den=median(Catch_KG), count=length(Catch_KG),
             lwr=quantile(Catch_KG, probs=.1),
             upr=quantile(Catch_KG, probs=.9))
xx <- merge(tmp, Inputs$loc_x, by.x='knot_i', by.y='knot_x')
xx <- subset(xx, count>5)
ggplot(xx, aes(Year, y=den, ymin=lwr, ymax=upr, color=Gear)) +
  geom_linerange() + geom_point()+ facet_wrap('knot_i') + scale_y_log10() +
  geom_vline(xintercept=2011)



## Try profiling all of the fixed effects. Maybe something with pop up?
source("startup.R")
## Test combined spatial model
n_x <- 50
model <- 'combined'; space <- 'ST'
savedir <- paste0(getwd(), '/fit_', model, "_", space,  "_", n_x)
set.seed(111) ## seed 111 works for ST; 112 crashes out
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=FALSE, loopnum=1,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=0, control=list(iter.max=300, trace=1))
mle <- Opt$par
get.profile <- function(i, mle){
  set.seed(111) ## seed 111 works for ST; 112 crashes out
  n_x <<- 50
  model <<- 'combined'; space <<- 'ST'
  combinedoff <<- FALSE
  savedir <<- paste0(getwd(), '/testing/test_prof_', i)
  source('startup.R')
  source("prepare_inputs.R")
  Obj$par <- mle
  err <- tryCatch(prof <- tmbprofile(obj=Obj, name=i, lower=TmbList$Lower[i], upper=TmbList$Upper[i]),
                  error=function(e) NULL)
  if(is.null(err)){
    ## Crashed out
    out <- list(crashed=TRUE, par=names(Obj$par)[i], prof=NULL)
  } else {
    out <- list(crashed=FALSE, par=names(Obj$par)[i], prof=prof)
  }
  return(out)
}
library(snowfall)
cores <- 18
pars <- 1:length(Obj$par)
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='prof_progress.txt')
sfExport('get.profile')
out.parallel <- sfLapply(pars, function(i) get.profile(i, Obj))
saveRDS('testing/prof.parallel.RDS', object=out.parallel)
sfStop()
par(mfrow=c(5,4))
for(i in 1:length(pars))
  err <- tryCatch(plot.tmbprofile(out.parallel[[i]]$prof),
                  error=function(e) NULL)


