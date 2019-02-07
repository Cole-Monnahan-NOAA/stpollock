## script to dig deeper into convergence issues by looking deeper at betas ñ

source("startup.R")
## L_beta1_z    L_beta1_z    L_beta1_z Beta_mean1_c Beta_mean1_c Beta_mean1_c
##      0.21         0.30         0.50         0.69        -0.29        -0.48
## L_beta2_z    L_beta2_z    L_beta2_z Beta_mean2_c Beta_mean2_c Beta_mean2_c
##      0.59         0.25         0.83         6.41         5.47         6.80
## logSigmaM    logSigmaM    logSigmaM
##      1.02         0.97         1.02


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
