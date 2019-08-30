## Explore the temporal cnofiguration for the AT data set using MLE models
## (For now at least).

setwd('..')
source("startup.R")
dir.create('sensitivities/temporal')

## Repeat with just the ATS
indices <- list(); k <- 1
for(model in c('ats', 'bts')){
  for(type in c('AR1', 'RW')){
    for(space in c('NS',"S", "ST")){
      for(extrayears in c(TRUE,FALSE)){
        control <- list(model=model, n_x=100)
        control$n_eps1 <- control$n_eps2 <- ifelse(space=='ST', 1,0)
        control$n_omega1 <- control$n_omega2 <- ifelse(space=='NS', 0,1)
        control$temporal <- ifelse(type=='AR1',4 ,2)
        control$replicateyears <- extrayears
        savedir <- paste0(getwd(), '/sensitivities/temporal/senfit_',model, "_", type,"_",
                          space)
        if(extrayears) savedir <- paste0(savedir, "_extrayears")
        source("prepare_inputs.R")
        Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                        upper=TmbList$Upper,   savedir=savedir,
                        newtonsteps=1, control=list(trace=1))
        results <- process.results(Opt, Obj, Inputs, model, space, savedir)
        indices[[k]] <- cbind(extrayears=extrayears, type=type, results$Index)
        k <- k+1
        plot.vastfit(results, plotmaps=TRUE)
      }
    }
  }
}

out <- do.call(rbind, indices)
out$space <- factor(out$space, levels=c('NS', 'S', 'ST'),
       labels=c('No space', 'Spatial', 'Spatiotemporal'))
out$extrayears <- factor(out$extrayears, levels=c(TRUE,FALSE),
                         labels=c("Extra years", "Not extra years"))
saveRDS(out, file='results/sensitivity_temporal.RDS')

g <- ggplot(subset(out, model=='ats'), aes(year, y=est, ymin=lwr, ymax=upr, color=type, fill=type)) +
  geom_ribbon(alpha=.5) + geom_line() + facet_grid(space~extrayears) + theme_bw()
ggsave('plots/sensitivity_temporal_ats.png', g, width=7, height=5)
g <- ggplot(subset(out, model=='bts'), aes(year, y=est, ymin=lwr, ymax=upr, color=type, fill=type)) +
  geom_ribbon(alpha=.5) + geom_line() + facet_grid(space~extrayears) + theme_bw()
ggsave('plots/sensitivity_temporal_bts.png', g, width=7, height=5)
