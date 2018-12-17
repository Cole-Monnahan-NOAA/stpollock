## File to run the fits to the real data
source('startup.R')
## Models are (combined, bts only, ats only) x (no space, space, spatiotemporal)


n_x <- 400 # number of knots
for(m in 1:3){
for(s in 2:3){
model <- c('ats', 'bts', 'combined')[m]
space <- c('NS', 'S', 'ST')[s]
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower,
                upper=TmbList$Upper,  savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## TMBhelper::Check_Identifiable(Obj)
results <- process.results(Opt, Obj, model, space, savedir)
plot.vastfit(results)
}
}


## Process them
xx <- list.dirs(, full.names=FALSE, recursive=FALSE)
dirs <- xx[grep('VAST_', xx)]

indices <- ldply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    return(Save$Index)
  } else {
    return(NULL)
  }
})
ggplot(indices, aes(year, est, group=strata, color=strata,  fill=strata)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.33) +
  geom_line() + geom_point()+
  facet_grid(model~space)

ests <- ldply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    m <- Save$Index$model[1]; s <- Save$Index$space[1]
    est <- data.frame(model=m, space=s, par=names(Save$Opt$SD$par.fixed),
                      est=Save$Opt$SD$par.fixed,
                      se=sqrt(diag(Save$Opt$SD$cov.fixed)), stringsAsFactors=FALSE)
    ## Add a number after each par to make them unique for plotting
    tmp <- as.vector(do.call(c, sapply(unique(est$par), function(x) {
      y <- est$par[est$par==x]
      1:length(y)})))
    est$par2 <- paste(est$par, tmp, sep='_')
    est$parnum <- tmp+switch(as.character(m), ats=1/3, bts=2/3, combined=0) +
      ifelse(s=='S', 0,.5)
    return(est)
  } else {
    return(NULL)
  }
})
ggplot(ests, aes(parnum, est, color=model, shape=space)) +
  geom_point(size=3) + facet_wrap('par', scales='free') +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
