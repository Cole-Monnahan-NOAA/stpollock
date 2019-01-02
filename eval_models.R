## Read in results and evalulate the models




## Pull out the results for different quantities
library(plyr)
xx <- list.dirs(, full.names=FALSE, recursive=FALSE)
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

fields <- llply(dirs, function(x) {
  ff <- file.path(x, 'Save.RData')
  ## for now just getting the combined fits
  if(length(grep('combined', x=x))>0 & file.exists(ff)){
    load(ff)
    kappa <- exp(Save$Opt$SD$par.fixed['logkappa1'])
    return(data.frame(model=Save$Index$model[1], space=Save$Index$space[1],
                      omegainput=Save$Report$Omegainput1_sf,
                      omega=Save$Report$Omega1_sc))
  } else {
    return(NULL)
  }
})
fields <- do.call(rbind.fill, fields)
fields <- cbind(fields, loc)
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
  scale_colour_gradientn(colours = Col(15))
ggsave('plots/map_omegas.png', width=9, height=6, units='in')

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


library(ggplot2)
g <- ggplot(indices, aes(year, est, group=model, color=model,  fill=model)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.33) +
  geom_line() + geom_point()+
  facet_grid(strata~space)
ggsave('plots/initial_fits.png', g, width=7, height=5)
## g <- ggplot(subset(ests, par=='beta1_ct'), aes(parnum, est, color=model, shape=space)) +
##   geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
##   geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
## ggsave('plots/initial_fits_beta1.png', g, width=7, height=5)
g <- ggplot(subset(ests, par!='beta1_ct'), aes(parnum, est, color=model, shape=space)) +
  geom_point(size=3) + facet_wrap('par', scales='free', ncol=2) +
  geom_linerange(aes(ymin=est-2*se, ymax=est+2*se), lwd=1.5)
ggsave('plots/initial_fits_pars.png', g, width=7, height=5)


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
