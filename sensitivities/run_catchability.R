
### Sensitivities to the structure of catchability coefficients lambda. The
### base case is a single value for lambda1. Sensitivities are (1)
### time-varying lambda (with tight prior) and (2) both lambda1 and lambda2
### constant
setwd(here())
chains <- 6
options(mc.cores = chains)
source('startup.R')
td <- 15
ad <- .9
iter <- 800
warmup <- 400

## The flag "fixlambda" controls the catchability
## configuration. Kind of a bad way to set it up, see
## prepare.inputs for how it actually works.

## Case 1: base case from paper, constant in p2
control <- list(model='combined', n_x=15, fixlambda=1,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/catchability_lambda2')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Case 2: time-varying catchability in P2.
control <- list(model='combined', n_x=15, fixlambda=-2,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/catchability_tvlambda2')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)

## Case 3: Constant in p1
control <- list(model='combined', n_x=15, fixlambda=2,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                make_plots=FALSE)
savedir <- paste0(getwd(), '/sensitivities/catchability_lambda1')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=iter, open_progress=FALSE, warmup=warmup,
               init=prior.fn, thin=1,
               control=list(max_treedepth=td, adapt_delta=ad))
saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
plot.mcmc(Obj, savedir, fit)


x1 <- readRDS('sensitivities/catchability_tvlambda2/results.mcmc.RDS')$index.strata
x2 <- readRDS('sensitivities/catchability_lambda2/results.mcmc.RDS')$index.strata
x3 <- readRDS('sensitivities/catchability_lambda2/results.mcmc.RDS')$index.strata
out <- rbind(cbind(x1, lambda='tv_lambda2'), cbind(x2, lambda='lambda2'),
             cbind(x3, lambda='lambda1'))
g <- ggplot(out, aes(year, est, ymin=lwr, ymax=upr, color=lambda, fill=lambda)) +
 geom_line(lwd=2) +  geom_ribbon(alpha=1/3) +
  facet_wrap('stratum', ncol=1) + theme_bw()
ggsave('plots/sensitivities_catchability.png', g, width=7, height=7)
