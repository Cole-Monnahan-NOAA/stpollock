## script to dig deeper into convergence issues by looking deeper at betas ñ

source("startup.R")
## L_beta1_z    L_beta1_z    L_beta1_z Beta_mean1_c Beta_mean1_c Beta_mean1_c
##      0.21         0.30         0.50         0.69        -0.29        -0.48
## L_beta2_z    L_beta2_z    L_beta2_z Beta_mean2_c Beta_mean2_c Beta_mean2_c
##      0.59         0.25         0.83         6.41         5.47         6.80
## logSigmaM    logSigmaM    logSigmaM
##      1.02         0.97         1.02


## Set seed and find a bad model without the ocmbination of strata and send
## it to Jim

## Try random seeds with combined model to test stability

## Test MCMC use bounds on L's given "close to MLE" crash point and run
## everything.


## Test combined spatial model
n_x <- 50
model <- 'combined'; space <- 'ST'
savedir <- paste0(getwd(), '/fit_', model, "_", space,  "_", n_x)
set.seed(112) ## seed 111 works for ST; 112 crashes out
options(warn=0)
source("prepare_inputs.R")
options(warn=2) # stop on a warning
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=TRUE, loopnum=3,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=0, control=list(iter.max=300, trace=1))
options(warn=0)
## Check where problem occurred
all <- Obj$env$last.par
fixed <- all[-Obj$env$random]
Obj$fn(fixed)
Obj$gr(fixed)
Obj$report(all)$jnll
hes <- Obj$env$spHess(par=all, random=TRUE)
evs <- eigen(hes)
plot(log10(evs$values+1e-4))

plot(Report$LogProb1_i)
plot(Report$LogProb2_i)
plot(log(Report$R2_i))
plot(log(Report$D_gcy))
plot(Report$Epsilon1_gct)
matplot(Report$Omegainput1_sf)
matplot(Report$Omegainput2_sf)
matplot(Report$beta1_tc)
matplot(Report$beta2_tc)

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
## This one is different
Map$beta2_ft <- factor(as.numeric(TmbList0$Map$beta2_ft) * NA)
Params$beta2_ft <- matrix(fixed[grep('beta2_ft', x=names(fixed))][TmbList0$Map$beta2_ft],nrow=3)
Map$logSigmaM <- factor(Params$logSigmaM*NA)
TmbList0$Parameters$logSigmaM
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

for(i in 1:length(Map)){
  print(names(Map)[i])
  print(length(Params[[names(Map)[i]]]))
  print(length(Map[[i]]))
#  print(length(TmbList0$Map[[i]]))
}

TmbList <- make_model(TmbData=TmbData, RunDir=savedir,
                      Version=Version,  RhoConfig=RhoConfig,
                      loc_x=Spatial_List$loc_x, Method=Method,
                      Param=Params, TmbDir='models',
                      Random=NULL, Map=Map)
Obj2  <-  TmbList[["Obj"]]
Obj2$env$beSilent()
Obj2$fn()
grs <- as.numeric(Obj2$gr())
plot(grs[1:50])
abline(v=36)

## Try optimizing just out of curiosity
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=FALSE,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=0, control=list(trace=20))




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
               iter=1000, open_progress=FALSE,
               init=inits, control=list(max_treedepth=7))
saveRDS(object = fit, file='fit.RDS')
fit <- readRDS('fit.RDS')
launch_shinystan(fit)



## Try to find a standard model that doesn't work to show Jim.
run.iteration <- function(seed){
  set.seed(seed)
  n_x <<- 50
  model <<- 'combined'; space <<- 'ST'
  combinedoff <<- FALSE
  savedir <<- paste0(getwd(), '/test_std_', seed)
  source('startup.R')
  source("prepare_inputs.R")
  options(warn=2) # stop immediately on NaN warning to save time
  err <- tryCatch(Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=FALSE,
                                  loopnum=3,
                                  upper=TmbList$Upper,  savedir=savedir,
                                  newtonsteps=0, control=list(iter.max=120, trace=1)),
                  error=function(e) NULL)
  if(is.null(err)){
    return('failed')
  } else {
    return(Opt[c( 'max_gradient', 'objective', 'par')])
    }
}

library(snowfall)
cores <- 10
chains <- cores*30
sfStop()
snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='convergence_progress.txt')
snowfall::sfExportAll()
sfLibrary(VAST)
sfLibrary(TMB)
out.parallel <-
  sfLapply(1:chains, function(i)
  ##lapply(1:chains, function(i)
  run.iteration(i))
sfStop()


which(out.parallel=='failed')
mean(out.parallel=='failed')
plot(sapply(out.parallel[out.parallel!='failed'], function(x) {x$max_gradient}))
plot(sapply(out.parallel[out.parallel!='failed'], function(x) {x$objective}))
write.csv(file='convergence.csv', x=t(sapply(out.parallel[out.parallel!='failed'],
            function(x) rbind(x$objective, x$max_gradient))))


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





## ## Loop through optimization one step at a time and save gradients and
## ## parameter vectors.
## pars <- grads <- nlls <- list()
## pars[[1]] <- Obj$par
## grads[[1]] <-  Obj$gr(Obj$par)
## nlls[[1]] <- Obj$fn(Obj$par)
## for(i in 2:100){
##   tmp <- optim(par=pars[[i-1]], fn=Obj$fn,
##                 gr=Obj$gr, lower=TmbList$Lower,
##                 upper=TmbList$Upper, method='L-BFGS-B',
##                 control=list(maxit=0, trace=1))
##   ## tmp <- nlminb(start=pars[[i-1]], objective=Obj$fn,
##   ##               gradient=Obj$gr, lower=TmbList$Lower,
##   ##               upper=TmbList$Upper,  savedir=savedir,
##   ##               control=list(iter.max=1, trace=0))
##   pars[[i]] <- tmp$par
##   grads[[i]] <- Obj$gr(tmp$par)
##   nlls[[i]] <- Obj$fn(tmp$par)
##   if(tmp$convergence==0) break
##   print(paste(i, round(max(grads[[i]]),4)))
## }
## pars <- as.data.frame(do.call(rbind,pars))
## grads <- as.data.frame(do.call(rbind,grads))
## nlls <- do.call(c, nlls)
## maxgrads <- apply(grads, 1, function(x) max(abs(x)))
## names(Obj$par)[apply(grads, 1, which.max)]
## names(pars) <- names(grads) <-
##   paste0(1:length(Obj$par),"_", names(Obj$par))
## grads$iter <- pars$iter <- 1:nrow(grads)
## pars.long <- melt(pars[-1,], id.vars='iter')
## grads.long <- melt(grads[-1,], id.vars='iter')
## ## Make some quick plots of these
## g <- ggplot(pars.long, aes(iter, value, group=variable, color=variable)) + geom_line()
## ggsave('testing/pars_NS_BTS.png', g, width=7, height=5)
## g <- ggplot(pars.long, aes(iter, value, group=variable)) +
##   geom_line() + facet_wrap('variable', scales='free_y')
## ggsave('testing/pars_faceted_NS_BTS.png', g, width=10, height=5)
## g <- ggplot(grads.long, aes(iter, value, group=variable, color=variable)) + geom_line()
## ggsave('testing/grads_NS_BTS.png', g, width=7, height=5)
## g <- ggplot(grads.long, aes(iter, sign(value)*log10(abs(value)), group=variable)) +
##   geom_line() + facet_wrap('variable')
## ggsave('testing/log_grads_NS_BTS.png', g, width=7, height=5)
## plot(nlls)
## plot(log10(maxgrads))
## prof <- tmbprofile(Obj, name=12)
## plot(prof)
## Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=FALSE,
##                 upper=TmbList$Upper,  savedir=savedir,
##                 newtonsteps=0, control=list(trace=1))
## Opt <- nlminb(start=Obj$par, objective=Obj$fn,
##               gradient=Obj$gr, lower=TmbList$Lower,
##               upper=TmbList$Upper,
##               control=list(iter.max=150, trace=0))

## library(tracer)
## library(nloptr)
## eval_f <- function(x) {
##   return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
## }
## ## Gradient of Rosenbrock Banana function
## eval_grad_f <- function(x) {
##   return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
##             200 * (x[2] - x[1] * x[1]) ) )
## }
## x0 <- c( -1.2, 1 )
## res <- nloptr( x0=x0,
##               eval_f=eval_f,
##               eval_grad_f=eval_grad_f,
##               opts=opts)
## res <- nloptr_tr( x0=x0,
##               eval_f=eval_f,
##               eval_grad_f=eval_grad_f,
##               opts=opts)
## print(res)
## ## tracer(res)
## opts <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8, print_level=2,"check_derivatives" = TRUE,
##              "check_derivatives_print" = "all")
## fn <- function(x) Obj$fn(x)
## gr <- function(x) Obj$gr(x)
## res <- nloptr(x0=Obj$par, eval_f=fn, eval_grad_f=gr, opts=opts,
##               lb=TmbList$Lower, ub=TmbList$Upper)
## print(res)









## Look at data by knot
tmp <- ddply(Data_Geostat, .(knot_i, Gear, Year), summarize,
            den=median(Catch_KG), count=length(Catch_KG))
xx <- subset(merge(tmp, Inputs$loc_x, by.x='knot_i', by.y='knot_x'),
             Gear !='Acoustic_16-surface' & ! Year %in% c(2011, 2013, 2015,
            2017))
Col  <-  colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
ggplot(xx, aes(E_km, N_km, col=log(den), size=sqrt(count), shape=den==0)) + geom_point() + facet_grid(Gear~Year)+
  scale_colour_gradientn(colours = Col(15)) + theme_bw()
x2 <- dcast(xx, knot_i+Year+E_km+N_km~Gear, value.var='den')
x3 <- dcast(xx, knot_i+Year+E_km+N_km~Gear, value.var='count')
names(x2)[5:6] <- names(x3)[5:6] <- c("BTS", "ATS")
x2$count.all <- x3$BTS + x3$ATS
ggplot(x2, aes(E_km, N_km, col=BTS>ATS, size=sqrt(count.all))) + geom_point() + facet_wrap(.~Year)+
  theme_bw()
ggplot(x2, aes(log(BTS),log(ATS))) + geom_point() + geom_abline(slope=1) + facet_wrap('Year')


## raw standardization
tmp <- subset(Data_Geostat, Catch_KG>0)
xx <- ddply(tmp, .(Gear, Year), summarize,
            catch=median(Catch_KG), std=sd(Catch_KG))
ggplot(xx, aes(Year, catch, group=Gear, color=Gear)) + geom_line() + geom_point()





### Old tests
## Test combined spatial model
n_x <- 50
model <- 'combined'; space <- 'S'
savedir <- paste0(getwd(), '/fit_', model, "_", space,  "_", n_x)
source("prepare_inputs.R")
## Look at report before inner optimization
r1 <- Obj$report()
pars1 <- matrix(Obj$env$last.par[grep('beta2_ft', names(Obj$env$last.par))],
               ncol=3, byrow=TRUE)
## Repeat after inner optimization
Obj$fn()
r2 <- Obj$report()
pars2 <- matrix(Obj$env$last.par[grep('beta2_ft', names(Obj$env$last.par))],
               ncol=3, byrow=TRUE)
Obj$env$beSilent()
## Try setting lower bound on the means
TmbList$Lower[grep('Beta_mean2_c', names(TmbList$Lower))]  <- 3.5
TmbList$Upper[grep('Beta_mean2_c', names(TmbList$Upper))]  <- 6.5
TmbList$Lower[grep('L_beta2_z', names(TmbList$Lower))]  <- 0
TmbList$Upper[grep('L_beta2_z', names(TmbList$Upper))]  <- 1
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, getsd=TRUE, loopnum=5,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## Now after outer optimization
r3 <- Obj$report()
pars3 <- matrix(Obj$env$last.par[grep('beta2_ft', names(Obj$env$last.par))],
               ncol=3, byrow=TRUE)


png(file.path(savedir, 'beta_tests.png'), width=7, height=6, res=400, units='in')
par(mfrow=c(3,3), mar=c(3,3,1.5,.5), oma=c(0,0,0,0), mgp=c(1.5,.5,0))
ylim1 <- log(range(c(r1$Index_cyl[,,1], r2$Index_cyl[,,1],r3$Index_cyl[,,1])))
ylim2 <- c(-4,4)
ylim3 <- range(c(r1$beta2_tc,r2$beta2_tc,r3$beta2_tc))
matplot(x=years, pars1, ylab='beta2_ft', type='l', ylim=ylim2,
        main='Pre-inner optim')
matplot(x=years, pars2, ylab='beta2_ft', type='l', ylim=ylim2, main='Post inner optim')
matplot(x=years, pars3, ylab='beta2_ft', type='l', ylim=ylim2, main='Post outer optim')
matplot(x=years, r1$beta2_tc, type='l', ylab='beta2_ct', ylim=ylim3)
matplot(x=years, r2$beta2_tc, type='l', ylab='beta2_ct', ylim=ylim3)
matplot(x=years, r3$beta2_tc, type='l', ylab='beta2_ct', ylim=ylim3)
matplot(x=years, y=log(t(r1$Index_cyl[,,1])), ylab='Log Index', type='l',
        ylim=ylim1)
matplot(x=years, y=log(t(r2$Index_cyl[,,1])), ylab='Log Index', type='l', ylim=ylim1)
matplot(x=years, y=log(t(r3$Index_cyl[,,1])), ylab='Log Index', type='l', ylim=ylim1)
dev.off()


## Look at P2 and R2
yy <- data.frame(year=Data_Geostat$Year, P2.1=r3$P2_iz[,1],
                 P2.2=r3$P2_iz[,2], gear=Data_Geostat$Gear)
yy <- melt(yy, id.vars=c('year', 'gear'))
g <- ggplot(yy, aes(y=value, x=factor(year))) + geom_violin() +
  facet_grid(gear~variable)
ggsave(file.path(savedir, 'P2_tests.png'), g, width=7, height=5)
xx <- data.frame(year=Data_Geostat$Year, R2=r3$R2_i, gear=Data_Geostat$Gear)
g <- ggplot(xx, aes(y=R2, x=factor(year))) + geom_violin() +
 facet_wrap('gear', ncol=1) + scale_y_log10()
ggsave(file.path(savedir, 'R2_tests.png'), g, width=7, height=5)

## Predicted vs observed for both models. Note that r4 is the copy of r3
## with no bounds on it. So you have to rerun that manually to recreate
## this
df <- data.frame(obs=Data_Geostat$Catch_KG, initial=r1$R2_i, unconstrained=r3$R2_i,
                 constrained=r4$R2_i, gear=Data_Geostat$Gear)
df <- subset(df, obs>0)
df <- melt(df, id.vars=c('gear','obs'), value.name='predicted')
g <- ggplot(df, aes(log(obs), log(predicted), color=gear)) + facet_grid(gear~variable) + geom_point(alpha=.5) +
  geom_abline(slope=1, intercept=0)
ggsave(file.path(savedir, 'obs_vs_pred_tests.png'), g, width=7, height=5)

## profile on beta meancz
prof1 <- tmbprofile(Obj, 24, parm.range=c(0,7), ystep=.001, h=1e-01, ytol=10)
prof1 <- tmbprofile(Obj, 24)
plot(prof1)
## Profile on variance of beta2[1]
prof2 <- tmbprofile(Obj, 20)
plot(prof)


### Some MCMC tests
library(tmbstan)
library(shinystan)
options(mc.cores = 7)
lwr <- TmbList$Lower
## L_betas are variances so bound below
lwr[grep('L_beta', names(lwr))] <- 0
mcmc <- tmbstan(obj=Obj, iter=1000, chains=7,
                init='par', upper= TmbList$Upper,
                lower=lwr,
                 control=list(max_treedepth=12, adapt_delta=.8),
                open_progress=FALSE)
saveRDS(mcmc, file='mcmc.RDS')
launch_shinystan(mcmc)

pars <- names(mcmc)[-grep('Omega', x = names(mcmc))]
pars <- names(mcmc)[grep(x=names(mcmc), 'beta1_ft')][(1:12)*3-1]
png('pairs.png', width=12, height=8, res=500, units='in')
pairs(mcmc, pars=pars)
dev.off()
