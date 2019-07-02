### Fit some resolution tests on the two independent model since they're
### stable
library(snowfall)
sfInit(parallel=TRUE, cpus=8)
inputs <- expand.grid(n_x=seq(210, 400, by=10), fs=c(FALSE, TRUE)[2], model=c('ats', 'bts'),
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
    Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=5, getsd=TRUE,
                    upper=TmbList$Upper,   savedir=savedir,
                    newtonsteps=0, control=list(trace=10)),
    error=function(e) 'xx')
  if(result=='xx'){
    cat('Failed at', ii, '\n')
    return(NULL)
  } else {
    results <- process.results(Opt, Obj, Inputs, inputs$model[ii], space, savedir)
    results$n_x <- n_x; results$finescale <- fs
    saveRDS(results, file=paste0(savedir, '/Save.RDS'))
    return(NULL)
  }
})




## After running everything, read the results back in and process them for
## plotting.
setwd('resolution_tests')
dirs <- list.dirs(getwd())[-1]
results.list <- lapply(dirs, function(x) {
  setwd(x)
  if(file.exists('Save.RDS'))
    Save <- readRDS('Save.RDS') else Save  <- NULL
  setwd('..')
  return(Save)
})
setwd('..')


library(ggplot2)
library(tidyr)
library(plyr)
pars <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
  data.frame(x$est, par2=paste0(1:nrow(x$est),'_', x$est$par),
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
g <- ggplot(pars, aes(n_x, y=est, color=model, lty=finescale)) +
#  geom_linerange(aes(ymin=lwr, ymax=upr), alpha=.5) +
  facet_wrap('par', scales='free_y') + geom_line(lwd=1) + theme_bw() +
  geom_vline(xintercept=50)
ggsave('plots/resolution_pars.png', g, width=12, height=5)

indices <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
  data.frame(x$Index.strata, n_x=ifelse(x$finescale, x$n_x, x$n_x+5),
             finescale=x$finescale, stringsAsFactors=FALSE)))
g <- ggplot(subset(indices, model=='ats'), aes(n_x, y=est, color=finescale)) +
  geom_linerange(aes(ymin=lwr, ymax=upr), alpha=.5) +
  facet_wrap('year', scales='free_y') + geom_line(lwd=1) + theme_bw() +
  geom_vline(xintercept=50) + ylab('log index')
ggsave('plots/resolution_indices_ats.png', g, width=12, height=5)
g <- ggplot(subset(indices, model=='bts'), aes(n_x, y=est, color=finescale)) +
  geom_linerange(aes(ymin=lwr, ymax=upr), alpha=.5) +
  facet_wrap('year', scales='free_y') + geom_line(lwd=1) + theme_bw() +
  geom_vline(xintercept=50) + ylab('log index')
ggsave('plots/resolution_indices_bts.png', g, width=12, height=5)


opts <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
          data.frame(nll=x$Opt$objective,
                     runtime=as.numeric(x$Opt$time_for_MLE, 'mins'),
                     sdtime=as.numeric(x$Opt$time_for_sdreport, 'mins'),
             iterations=x$Opt$iterations, evaluations=x$Opt$evaluations[1],
             maxgrad=x$Opt$max_gradient, AIC=x$Opt$AIC,
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
opts.long <- gather(opts, variable, value, -n_x, -finescale, -model)
g <- ggplot(opts.long, aes(n_x, y=value, color=model, lty=finescale)) +
  facet_wrap('variable', scales='free_y') + geom_line() + theme_bw() +
  geom_vline(xintercept=50)
ggsave('plots/resolution_opts.png', width=9, height=5)

nlls <- do.call(rbind, lapply(results.list, function(x) if(!is.null(x))
  data.frame(comp=c('omega1', 'eps1', 'omega2', 'eps2', 'vessel1',
                    'vessel2', '6', '7', 'betapen1', 'betapen2', 'data1', 'data2',
                    'delta', 'Bstruct', 'SVcov1', 'SVcov2'),
             nll=x$Report$jnll_comp,
             model=x$model, n_x=ifelse(x$model=='bts', x$n_x, x$n_x+3),
             finescale=x$finescale, stringsAsFactors=FALSE)))
nlls <- ddply(subset(nlls, nll>0), .(model, comp), mutate,
              nll0=nll-nll[1])
g <- ggplot(nlls, aes(n_x, y=nll0, color=model, lty=finescale)) +
  facet_wrap('comp', scales='free_y') + geom_line() + theme_bw()
ggsave('plots/resolution_nllcomp.png', g, width=9, height=5)
