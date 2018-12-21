## Read in results and evalulate the models




## Pull out the results for different quantities
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

fields <- llply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  ## for now just getting the combined fits
  if(length(grep('combined', x=x))>0 & file.exists(ff)){
    load(ff)
    return(data.frame(model=Save$Index$model[1], space=Save$Index$space[1],
                      omegainput=Save$Report$Omegainput1_sf,
                      omega=Save$Report$Omega1_sc))
  } else {
    return(NULL)
  }
})
fields <- do.call(rbind.fill, fields)

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


g <- ggplot(indices, aes(year, est, group=strata, color=strata,  fill=strata)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.33) +
  geom_line() + geom_point()+
  facet_grid(model~space)
ggsave('plots/initial_fits.png', g, width=7, height=5)
g <- ggplot(subset(ests, par=='beta1_ct'), aes(parnum, est, color=model, shape=space)) +
  geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
ggsave('plots/initial_fits_beta1.png', g, width=7, height=5)
g <- ggplot(subset(ests, par!='beta1_ct'), aes(parnum, est, color=model, shape=space)) +
  geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
ggsave('plots/initial_fits_pars.png', g, width=7, height=5)


## Get the correlation matrices out for the combined models
ff <- file.path('VAST_output_combined_ST', 'Save.RData')
load(ff)
L_vec <- Save$ParHat[names(Save$ParHat)=='L_omega1_z']
Loadings_pf <- matrix(NA, 3,3)
Count <- 1
for(f in 1:3){
  for(p in 1:3){
    if(p>=f){
      Loadings_pf[p,f] = L_vec[Count]
      Count <- Count + 1
    }else{
      Loadings_pf[p,f] = 0.0;
    }
  }
}

Cov <- t(Loadings_pf) %*% Loadings_pf
Cor <- round(cov2cor(Cov),3)

