
source("startup.R")

## Use the posterior medians for initial values of fixed effects
fit <- readRDS("mcmc_basecase_100/mcmcfit.RDS")
pars.fixed <- apply(as.data.frame(fit)[,-Obj$env$random], 2, median)
pars.all <- apply(as.data.frame(fit), 2, median)

## Base case for paper: combined model
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='combined',
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_100')
source("prepare_inputs.R")
Obj$env$last.par <- pars.all[-length(pars.all)]
Obj$par <- pars.fixed[-length(pars.fixed)]
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Repeat with just the BTS
control <- list(seed=121, beta2temporal=TRUE, n_x=100, model='bts',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_100_bts')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Repeat with just the ATS
control$model <- 'ats'
savedir <- paste0(getwd(), '/fit_basecase_100_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Base case for paper w/ finescale on
control <- list(seed=121, beta2temporal=TRUE, n_x=50,
                n_eps1="IID", n_eps2="IID", n_omega2="IID", n_omega1="IID",
                beta1temporal=TRUE, filteryears=FALSE, finescale=TRUE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=FALSE)
savedir <- paste0(getwd(), '/fit_basecase_finscale')
source("prepare_inputs.R")
Obj$par <- pars.fixed[-length(pars.fixed)]
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=6, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=0, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)


### Fit some resolution tests on the two independent model since they're
### stable
library(snowfall)
sfInit(parallel=TRUE, cpus=8)
inputs <- expand.grid(n_x=seq(210, 400, by=10), fs=c(FALSE, TRUE)[1], model=c('ats', 'bts'),
                      stringsAsFactors=FALSE)
sfExport('inputs')
trash <- sfLapply(1:nrow(inputs), function(ii){
  fs <<- inputs$fs[ii]
  n_x <<- inputs$n_x[ii]
  model <<- inputs$model[ii]
  source("startup.R")
  control <<- list(seed=121, beta2temporal=TRUE, n_x=n_x, finescale=fs,
                   n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                   ## n_eps1=0, n_eps2=0, n_omega2=1, n_omega1=0,
                   beta1temporal=TRUE, filteryears=FALSE,
                   kappaoff=0, temporal=2, fixlambda=12,
                   make_plots=FALSE, model=model)
  savedir <<- paste0(getwd(), '/resolution_tests/test_', n_x, '_', model)
  if(fs) savedir <<- paste0(savedir, '_finescale')
  source("prepare_inputs.R")
  result <- tryCatch(
    Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                    upper=TmbList$Upper,   savedir=savedir,
                    newtonsteps=0, control=list(trace=10)),
    error=function(e) 'xx')
  if(result=='xx'){
    return(NULL)
  } else {
    results <- process.results(Opt, Obj, Inputs, inputs$model[ii], space, savedir)
    results$n_x <- n_x; results$finescale <- fs
    saveRDS(results, file=paste0(savedir, '/Save.RDS'))
    return(NULL)
  }
})

setwd('resolution_tests')
dirs <- list.dirs(getwd())[-1]
results.list <- lapply(dirs, function(x) {
  setwd(x)
  print(x)
  Save <- readRDS('Save.RDS')
  setwd('..')
  return(Save)
})
setwd('..')


library(ggplot2)
library(tidyr)
library(plyr)
pars <- do.call(rbind, lapply(results.list, function(x)
  data.frame(x$est, par2=paste0(1:nrow(x$est),'_', x$est$par),
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
g <- ggplot(pars, aes(n_x, y=est, color=model)) + geom_linerange(aes(ymin=lwr, ymax=upr)) +
  facet_wrap('par', scales='free_y') + geom_line() + theme_bw()
ggsave('plots/resolution_pars.png', width=9, height=5)

opts <- do.call(rbind, lapply(results.list, function(x)
  data.frame(nll=x$Opt$objective, runtime=as.numeric(x$Opt$time_for_MLE),
             iterations=x$Opt$iterations, evaluations=x$Opt$evaluations[1],
             maxgrad=x$Opt$max_gradient, AIC=x$Opt$AIC,
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
opts.long <- gather(opts, variable, value, -n_x, -finescale, -model)
g <- ggplot(opts.long, aes(n_x, y=value, color=model)) +
  facet_wrap('variable', scales='free_y') + geom_line() + theme_bw()
ggsave('plots/resolution_opts.png', width=9, height=5)

nlls <- do.call(rbind, lapply(results.list, function(x)
  data.frame(comp=c('omega1', 'eps1', 'omega2', 'eps2', 'vessel1',
                    'vessel2', '6', '7', 'betapen1', 'betapen2', 'data1', 'data2',
                    'delta', 'Bstruct', 'SVcov1', 'SVcov2'),
             nll=x$Report$jnll_comp,
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
nlls <- ddply(subset(nlls, nll>0), .(model, comp), mutate,
              nll0=nll-nll[1])
g <- ggplot(nlls, aes(n_x, y=nll0, color=model)) +
  facet_wrap('comp', scales='free_y') + geom_line() + theme_bw()
ggsave('plots/resolution_nllcomp.png', width=9, height=5)
