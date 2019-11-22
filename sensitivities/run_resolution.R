## File to run the fits to the real data
chains <- 6
options(mc.cores = chains)
td <- 15
ad <- .8
iter <- 800
warmup <- 400
dir.create('sensitivities/resolution', showWarnings=FALSE)

## Base case for paper: combined
for(nx in c(50, 100, 200, 300)){
  control <- list(model='combined', n_x=nx,
                  n_eps1="IID", n_eps2="IID",
                  n_omega2="IID", n_omega1="IID",
                  make_plots=FALSE)
  savedir <- paste0(getwd(), '/sensitivities/resolution/senfit_', nx)
  source("prepare_inputs.R")
  fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                 iter=iter, open_progress=FALSE, warmup=warmup,
                 init=prior.fn, seed=12512,
                 control=list(max_treedepth=td, adapt_delta=ad))
  saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
  plot.mcmc(Obj, savedir, fit)
}

## Read in results and save them
dirs <- list.dirs('sensitivities/resolution', recursive=FALSE)
results.list <- lapply(dirs, function(x) {
  nx <- as.numeric(strsplit(x, '_')[[1]][2])
  res <- readRDS(file.path(x, 'results.mcmc.RDS'))
  index <- cbind(res$index.strata, knots=nx)
  return(index)
})
## Create a special x column to improve plotting look
results <- do.call(rbind, results.list) %>%
  mutate(knots2=as.numeric(factor(knots), levels=c(50,100,200,300))+ .15 * as.numeric(stratum))
saveRDS(results, 'results/resolution.RDS')

g <- ggplot(results, aes(knots2, est, ymin=lwr, ymax=upr, color=stratum)) +
  geom_linerange( lwd=1) +geom_line() + geom_point(size=2) +
  facet_wrap('year') + theme_bw() + ylab('log index') +
  scale_x_continuous(labels=sort(unique(results$knots)),
  breaks=1:4) + xlab('# of knots')
ggsave('plots/resolution_by_knots.png', g, width=9, height=5.5)

## stop("This is old")
## ### Fit some resolution tests on the two independent model since they're
## ### stable
## library(snowfall)
## sfInit(parallel=TRUE, cpus=8)
## inputs <- expand.grid(n_x=seq(210, 400, by=10), fs=c(FALSE, TRUE)[1], model=c('ats', 'bts'),
##                       stringsAsFactors=FALSE)
## sfExport('inputs')
## trash <- sfLapply(1:nrow(inputs), function(ii){
##   fs <<- inputs$fs[ii]
##   n_x <<- inputs$n_x[ii]
##   model <<- inputs$model[ii]
##   source("startup.R")
##   control <<- list(seed=121, beta2temporal=TRUE, n_x=n_x, finescale=fs,
##                    n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
##                    ## n_eps1=0, n_eps2=0, n_omega2=1, n_omega1=0,
##                    beta1temporal=TRUE, filteryears=FALSE,
##                    kappaoff=0, temporal=2, fixlambda=12,
##                    make_plots=FALSE, model=model)
##   savedir <<- paste0(getwd(), '/resolution_tests/test_', n_x, '_', model)
##   if(fs) savedir <<- paste0(savedir, '_finescale')
##   source("prepare_inputs.R")
##   result <- tryCatch(
##     Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
##                     upper=TmbList$Upper,   savedir=savedir,
##                     newtonsteps=0, control=list(trace=10)),
##     error=function(e) 'xx')
##   if(result=='xx'){
##     return(NULL)
##   } else {
##     results <- process.results(Opt, Obj, Inputs, inputs$model[ii], space, savedir)
##     results$n_x <- n_x; results$finescale <- fs
##     saveRDS(results, file=paste0(savedir, '/Save.RDS'))
##     return(NULL)
##   }
## })

## setwd('resolution_tests')
## dirs <- list.dirs(getwd())[-1]
## results.list <- lapply(dirs, function(x) {
##   setwd(x)
##   print(x)
##   if(file.exists('Save.RDS'))
##     Save <- readRDS('Save.RDS') else Save  <- NULL
##   setwd('..')
##   return(Save)
## })
## setwd('..')


## library(ggplot2)
## library(tidyr)
## library(plyr)
## pars <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
##   data.frame(x$est, par2=paste0(1:nrow(x$est),'_', x$est$par),
##              model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
##              finescale=x$finescale, stringsAsFactors=FALSE)))
## g <- ggplot(pars, aes(n_x, y=est, color=model)) + geom_linerange(aes(ymin=lwr, ymax=upr)) +
##   facet_wrap('par', scales='free_y') + geom_line() + theme_bw()
## ggsave('plots/resolution_pars.png', width=9, height=5)

## opts <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
##           data.frame(nll=x$Opt$objective,
##                      runtime=as.numeric(x$Opt$time_for_MLE, 'mins'),
##                      sdtime=as.numeric(x$Opt$time_for_sdreport, 'mins'),
##              iterations=x$Opt$iterations, evaluations=x$Opt$evaluations[1],
##              maxgrad=x$Opt$max_gradient, AIC=x$Opt$AIC,
##              model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
##              finescale=x$finescale, stringsAsFactors=FALSE)))
## opts.long <- gather(opts, variable, value, -n_x, -finescale, -model)
## g <- ggplot(opts.long, aes(n_x, y=value, color=model)) +
##   facet_wrap('variable', scales='free_y') + geom_line() + theme_bw()
## ggsave('plots/resolution_opts.png', width=9, height=5)

## nlls <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
##   data.frame(comp=c('omega1', 'eps1', 'omega2', 'eps2', 'vessel1',
##                     'vessel2', '6', '7', 'betapen1', 'betapen2', 'data1', 'data2',
##                     'delta', 'Bstruct', 'SVcov1', 'SVcov2'),
##              nll=x$Report$jnll_comp,
##              model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
##              finescale=x$finescale, stringsAsFactors=FALSE)))
## nlls <- ddply(subset(nlls, nll>0), .(model, comp), mutate,
##               nll0=nll-nll[1])
## g <- ggplot(nlls, aes(n_x, y=nll0, color=model)) +
##   facet_wrap('comp', scales='free_y') + geom_line() + theme_bw()
## ggsave('plots/resolution_nllcomp.png', width=9, height=5)


## Explore the raw indices and how it changes with resoultion
## for(nx in c(10, 25, 50, 75, 100, 200, 300, 400, 500, 750)){
##   control <- list(model='combined', n_x=nx,
##                   n_eps1="IID", n_eps2="IID",
##                   n_omega2="IID", n_omega1="IID",
##                   make_plots=TRUE)
##   savedir <- paste0(getwd(), '/sensitivities/resolution/senfit_', nx)
##   source("prepare_inputs.R")
##   ## fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
##   ##                iter=iter, open_progress=FALSE, warmup=warmup,
##   ##                init=prior.fn, seed=12512,
##   ##                control=list(max_treedepth=td, adapt_delta=ad))
##   ## saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
##   ## plot.mcmc(Obj, savedir, fit)
## }
## dirs <- list.dirs('sensitivities/resolution', recursive=FALSE)
## index.raw <- lapply(dirs, function(x) {
##   nx <- as.numeric(strsplit(x, '_')[[1]][2])
##   res <- read.csv(file.path(x, 'index.raw.csv'))
##   index.raw <- cbind(res, knots=nx)
##   return(index.raw)
## }) %>% do.call(rbind, .)
## g <- ggplot(index.raw, aes(year, log(value), group=knots, color=factor(knots))) +
##   geom_line() + geom_point(size=2) +
##   facet_grid(gear~type, scales='free_y') + theme_bw() + ylab('log index')
