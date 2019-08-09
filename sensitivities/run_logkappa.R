## A series of sensitivity analyses to run
chains <- 4
options(mc.cores = chains)
setwd('..')
source('startup.R')
td <- 15
ad <- .9
iter <- 800
warmup <- 200


### The effect of fixing logkappa. Run the models with half and double the
### spatial range used (50km and 200km)
for(model in c('bts', 'ats', 'combined')[-3]){
  for(kappascale in c(.5,1, 2)){
    control <- list(beta2temporal=TRUE, n_x=100,
                n_eps1=1, n_eps2=1, n_omega2=1, n_omega1=1,
                model=model, kappascale=kappascale,
                kappaoff=12)
    if(model=='combined')
      control[c('n_eps1', 'n_eps2', 'n_omega1', 'n_omega2')] <- "IID"
    savedir <- paste0(getwd(), '/sensitivities/mcmcfit_kappascale_', kappascale,'_', model)
    source("prepare_inputs.R")
    fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
                   iter=iter, open_progress=FALSE, warmup=warmup,
                   init=prior.fn, thin=1,
                   control=list(max_treedepth=td, adapt_delta=ad))
    saveRDS(object = fit, file=paste0(savedir,'/mcmcfit.RDS'))
    plot.mcmc(Obj, savedir, fit)
  }
}

### Make quick plots of these sensitivities
##The indices
results.list <- readRDS('results/sensitivity_aniso.RDS')
out <- do.call(rbind, lapply(results.list, function(x) data.frame(range=factor(100*x$kappascale), x$Index)))
g <- ggplot(out, aes(year, est, color=range)) + geom_line() +
  facet_wrap('model', ncol=1) + theme_bw() +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=range), alpha=.3) + ylab('log index')
ggsave('plots/sensitivity_aniso_indices.png', g, width=7, height=5)
## Aniso estimates
out <- do.call(rbind, lapply(results.list, function(x)
  data.frame(range=100*x$kappascale, x$est, model=x$model))) %>%
  filter(par=='ln_H_input') %>% cbind(par2=c('ln_H_1','ln_H_2'))
g <- ggplot(out, aes(range, est)) + geom_line() +
  facet_grid(par2~model, scales='free_y') + theme_bw()+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.3)
ggsave('plots/sensitivity_aniso_H_estimates.png', g, width=7, height=5)
## spatial variance estimates
out <- do.call(rbind, lapply(results.list, function(x)
  data.frame(range=100*x$kappascale, x$est, model=x$model))) %>%
  filter(grepl('L_eps|L_om', as.character(par)))
g <- ggplot(out, aes(range, est, group=par, color=par)) + geom_line() +
  facet_grid(.~model, scales='free_y') + theme_bw()
##  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.3)
ggsave('plots/sensitivity_aniso_L_estimates.png', g, width=7, height=5)
plot.aniso <- function(i, Range){
  Report <- results.list[[i]]$Report
  ## Stol this from plot_anisotropy so could plot all 4 together
  Eigen = eigen(Report$H)
  rss = function(V) sqrt(sum(V[1]^2 + V[2]^2))
  Pos_Major = Eigen$vectors[, 1] * Eigen$values[1] * Report$Range_raw1
  Pos_Minor = Eigen$vectors[, 2] * Eigen$values[2] * Report$Range_raw1
  Pres_Major = Eigen$vectors[, 1] * Eigen$values[1] * Report$Range_raw2
  Pres_Minor = Eigen$vectors[, 2] * Eigen$values[2] * Report$Range_raw2
  ## Range = 1.1 * c(-1, 1) * max(abs(cbind(Pos_Major, Pos_Minor, Pres_Major, Pres_Minor)))
  ##  print(Range)
  plot(1, type = "n", xlim = Range, ylim = c(Range[1], Range[2] * 1.2), xlab = "", ylab = "")
  shape::plotellipse(rx = rss(Pres_Major), ry = rss(Pres_Minor),
                     angle = -1 * (atan(Pres_Major[1]/Pres_Major[2])/(2 * pi) * 360 - 90), lcol = c("green", "black")[1],
                     lty = c("solid", "dotted")[1])
  shape::plotellipse(rx = rss(Pos_Major), ry = rss(Pos_Minor),
                     angle = -1 * (atan(Pos_Major[1]/Pos_Major[2])/(2 * pi)
                       * 360 - 90), lcol = "black", lty = "solid")
  mtext(paste0(results.list[[i]]$model, ', range=',
               100*results.list[[i]]$kappascale, ' km'), side=3, line=.5)
}
png('plots/sensitivity_aniso_ellipses.png', width=9, height=5, units='in', res=500)
par(mfrow=c(2,4), mar=c(3,3,2.5,.5))
for(i in 1:8) plot.aniso(i, Range=c(-6316, 6316))
## mtext(side = 1, outer = FALSE, line = 2, text = "Eastings (km.)")
## mtext(side = 2, outer = FALSE, line = 2, text = "Northings (km.)")
legend("top", legend = c("Encounter probability",
                         "Positive catch rates"), fill = c("green", "black"), bty = "n")
dev.off()
## Save a table of these
out <- do.call(rbind, lapply(results.list, function(x)
  data.frame(range=100*x$kappascale, x$est, model=x$model)))
write.csv(out, file='results/table.sensitivity.aniso.csv')
