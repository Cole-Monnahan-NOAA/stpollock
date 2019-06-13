## Read in results and evalulate the models
library(plyr)
library(ggplot2)
xx <- list.dirs(full.names=FALSE, recursive=FALSE)

## Pull out the results for the main fits to the model and plot key results
dirs <- xx[grep('fit_', xx)]
indices <- ldply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    return(Save$Index)
  } else {
    return(NULL)
  }
})
g <- ggplot(indices, aes(year, est, group=model, color=model,  fill=model)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.33) +
  geom_line() + geom_point()+
  facet_grid(strata~space)
ggsave('plots/initial_fits.png', g, width=7, height=5)

## Bias corrected test
dirs <- xx[grep('bias', xx)]
indices <- ldply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    y <- length(grep('cor', x))>0
    return(cbind(bias.corrected=y, Save$Index))
  } else {
    return(NULL)
  }
})
g <- ggplot(indices, aes(year, est, group=bias.corrected, color=bias.corrected,  fill=bias.corrected)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.33) +
  geom_line() + geom_point()+
  facet_grid(strata~space)
ggsave('plots/bias.corrected_test.png', g, width=7, height=5)


dirs2 <- dirs[-grep('NS', x=dirs)]
fields <- llply(dirs2, function(x) {
  ff <- file.path(x, 'Save.RData')
  ## for now just getting the combined fits
  if(length(grep('combined', x=x))>0 & file.exists(ff)){
    load(ff)
    kappa <- exp(Save$Opt$SD$par.fixed['logkappa1'])
    return(data.frame(model=Save$Index$model[1], space=Save$Index$space[1],
                      omegainput=Save$Report$Omegainput1_sf,
                      omega=Save$Report$Omega1_sc, E_km=Save$Inputs$loc$E_km,
                      N_km=Save$Inputs$loc$N_km))
  } else {
    return(NULL)
  }
})
fields <- do.call(rbind.fill, fields)
fields.long <- melt(fields, id.vars=c('model', 'space', 'E_km', 'N_km'),
                    factorsAsStrings=FALSE)
fields.long$strata <- paste0('strata_',unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
       function(x) x[2])))
fields.long$type <- unlist(lapply(strsplit(as.character(fields.long$variable), split='\\.'),
       function(x) x[1]))
fields.long$type <- factor(fields.long$type, levels=c('omegainput', 'omega'))
fields.long <- ddply(fields.long, .(type, space), mutate,
                     normalized=value/sd(value))
Col  <-  colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
g <- ggplot(fields.long, aes(E_km, N_km, col=normalized)) +
  geom_point(alpha=.5, size=1) +
  facet_grid(space+type~strata) +
  scale_colour_gradientn(colours = Col(15)) + theme_bw()
ggsave('plots/map_omegas.png', g, width=9, height=6, units='in')

ests <- ldply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    m <- Save$Index$model[1]; s <- Save$Index$space[1]
    est <- data.frame(model=m, space=s, par=names(Save$Opt$SD$par.fixed),
                      est=Save$Opt$SD$par.fixed,
                      se=sqrt(diag(Save$Opt$SD$cov.fixed)), stringsAsFactors=FALSE)
    ## Add a number after each par to make them unique for plotting
    tmp <- as.vector(do.call(c, lapply(unique(est$par), function(x) {
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
## g <- ggplot(subset(ests, par=='beta1_ct'), aes(parnum, est, color=model, shape=space)) +
##   geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
##   geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
## ggsave('plots/initial_fits_beta1.png', g, width=7, height=5)
g <- ggplot(ests, aes(parnum, est, color=model, shape=space)) +
  geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
g
ggsave('plots/initial_fits_pars.png', g, width=7, height=5)


##
load("fit_combined_NS/Save.RData")
beta1.NS <- data.frame(space='NS', Save$Report$beta1_tc)
load("fit_combined_S/Save.RData")
beta1.S <- data.frame(space='S', Save$Report$beta1_tc)
beta1 <- rbind(beta1.NS, beta1.S)
names(beta1) <- c('space', 'stratum1', 'stratum2', 'stratum3')
beta1$year.missing <-2007:2018 %in% c(2011, 2013, 2015, 2017)
beta1.long <- melt(beta1, id.vars=c('space', 'stratum3', 'year.missing'))
g <- ggplot(beta1.long, aes(y=stratum3, x=value, group=space, fill=space,
                            col=space, shape=year.missing)) +
  geom_abline(slope=1, intercept=0) +
  geom_point(size=1.5) +
  facet_grid(variable~space, scales='free') +
  stat_smooth(alpha=.25, method='lm') + theme_bw()
ggsave('plots/beta1_correlations.png', g, width=7, height=3)



### Resolution test
xx <- list.dirs(full.names=FALSE, recursive=FALSE)
dirs <- xx[grep('knots_', xx)]
ests <- llply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  if(file.exists(ff)){
    load(ff)
    m <- Save$Index$model[1]; s <- Save$Index$space[1]
    res <- as.numeric(strsplit(x, '_')[[1]][4])
    est <- data.frame(model=m, space=s, par=names(Save$Opt$SD$par.fixed),
                      est=Save$Opt$SD$par.fixed, res=res, AIC=as.numeric(Save$Opt$AIC),
                      se=sqrt(diag(Save$Opt$SD$cov.fixed)), stringsAsFactors=FALSE)
    ## Add a number after each par to make them unique for plotting
    tmp <- as.vector(do.call(c, lapply(unique(est$par), function(x) {
      y <- est$par[est$par==x]
      1:length(y)})))
    est$par2 <- paste(est$par, tmp, sep='_')
    est$parnum <- tmp+switch(as.character(m), ats=1/3, bts=2/3, combined=0) +
      ifelse(s=='S', 0,.5)
    index <- cbind(Save$Index, res=res)
    x <- list(est=est, index=index)
    return(x)
  } else {
    return(NULL)
  }
})
pars <- do.call(rbind, lapply(ests, function(x) x[[1]]))
indices <- do.call(rbind, lapply(ests, function(x) x[[2]]))
g <- ggplot(pars, aes(log2(res), est)) +
  geom_point() + geom_line() + facet_wrap('par2', scales='free_y', ncol=4) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se)) + xlab("log2(# knots)")
ggsave('plots/resolution_pars.png', g, width=7, height=5)
g <- ggplot(indices, aes(log2(res), est, group=strata, color=strata)) +
  geom_point() + geom_line() + facet_wrap('year', ncol=4) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se)) +
  xlab("log2(# knots)") + scale_y_log10()
ggsave('plots/resolution_indices.png', g, width=7, height=5)



## Get the correlation matrices out for the combined models
ff <- file.path('fit_combined_S', 'Save.RData')
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


## Mapdetails <- make_map_info(Region, NN_Extrap=Spatial_List$NN_Extrap,
##                             Extrapolation_List=Extrapolation_List)
## Plot_factors(Save$Report, Save$ParHat, Data=TmbData, SD=Save$Opt$SD,
##              mapdetails_list=Mapdetails)
## FishStatsUtils::PlotMap_Fn(
##              MappingDetails=mapdetails_list[["MappingDetails"]],
##              Mat=Mat_sf, PlotDF=mapdetails_list[["PlotDF"]],
##              MapSizeRatio=mapdetails_list[["MapSizeRatio"]],
##              Xlim=mapdetails_list[["Xlim"]],
##              Ylim=mapdetails_list[["Ylim"]],
##              Format='none',
##              ## FileName=paste0(plotdir,"Factor_maps--",Par_name),
##              Year_Set=paste0("Factor_",1:ncol(Mat_sf)),
##              Rotate=mapdetails_list[["Rotate"]],
##              zone=mapdetails_list[["Zone"]], mar=c(0,0,2,0),
##              oma=c(2.5,2.5,0,0), Cex=0.01, mfrow=Dim_factor, pch=20,
##              Legend=mapdetails_list[["Legend"]], plot_legend_fig=FALSE,
##              land_color=land_color)
