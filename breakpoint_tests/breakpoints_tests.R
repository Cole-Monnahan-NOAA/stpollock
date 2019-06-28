
## Look at how the ATS changes with different breakpoints for the lower and
## upper bounds (EFH)
setwd('..')

source("startup.R")

library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)

ats1 <- cbind(breakpoints='3_16', read.csv('data/ats_3_16.csv'))
ats2 <- cbind(breakpoints='.5_16', read.csv('data/ats_.5_16.csv'))
ats3 <- cbind(breakpoints='.5_20', read.csv('data/ats_.5_20.csv'))
ats4 <- cbind(breakpoints='.5_12', read.csv('data/ats_.5_12.csv'))
## ats2 <- ats2[,names(ats1)]
## ats3 <- ats3[,names(ats1)]
## ats4 <- ats4[,names(ats1)]

ats <- rbind(ats1, ats2, ats3, ats4) %>%
  select(-X, -surface, -ground, -dist, -strata1) %>%
  gather(stratum, value, -lon, -lat, -breakpoints, -year) %>%
  filter(value<10000)
g <- ggplot(ats, aes(log(value), fill=breakpoints)) +
  geom_histogram(position='identity', alpha=.5, bins=50) +
  facet_grid(year~stratum) + theme_bw()
ggsave('plots/breakpoints_histograms.png', g, width=7, height=5, dpi=500)
col  <-  colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
g <- ggplot(subset(ats, stratum=='strata2'), aes(lon ,lat, color=log(value))) +
  geom_point(cex=1) + facet_grid(year~breakpoints) +
  scale_color_gradientn(colours=col(15)) + theme_bw()
ggsave('plots/breakpoints_maps.png', g,  width=7, height=5, dpi=500)


controlc <- list(seed=121, beta2temporal=TRUE, n_x=200, model='combined',
                n_eps1='IID', n_eps2='IID', n_omega2='IID', n_omega1='IID',
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)

## Fit 3-16 (original)
control <- list(seed=121, beta2temporal=TRUE, n_x=200, model='ats',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_3_16_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_3_16_combined')
control <- controlc
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)


## Fit .5-16 !!!! Manually copy the ats.csv over !!!!!
control <- list(seed=121, beta2temporal=TRUE, n_x=200, model='ats',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_16_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_16_combined')
control <- controlc
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Fit .5-12 !!!! Manually copy the ats.csv over !!!!!
control <- list(seed=121, beta2temporal=TRUE, n_x=200, model='ats',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_12_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_12_combined')
control <- controlc
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)

## Fit .5-20 !!!! Manually copy the ats.csv over !!!!!
control <- list(seed=121, beta2temporal=TRUE, n_x=200, model='ats',
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                beta1temporal=TRUE, filteryears=FALSE, finescale=FALSE,
                kappaoff=12, temporal=2, fixlambda=2, make_plots=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_20_ats')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
savedir <- paste0(getwd(), '/breakpoint_tests/breakpoints_.5_20_combined')
control <- controlc
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=1))
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results, plotmaps=TRUE)
