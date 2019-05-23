source("startup.R")
model <- 'combined'

## In case any of them crash below
out1.1 <- out1.2 <- out2.1 <- out2.2 <- out3.1 <- out3.2 <-
  out4.1 <- out4.2 <- NULL
res1.1 <- res1.2 <- res2.1 <- res2.2 <- res3.1 <- res3.2 <-
  res4.1 <- res4.2 <- NULL


## Base case simplified model
control <- list(seed=121, beta2temporal=FALSE, beta1temporal=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0, n_x=100,
                combinedoff=FALSE, filteryears=TRUE, make_plots=TRUE,
                kappaoff=0, temporal=0, fixlambda=12)
savedir <- paste0(getwd(), '/test1_combined')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
index.raw <- read.csv(paste0(savedir, '/index.raw.csv'))
index.raw <- subset(index.raw, type=='Naive Spatial' & year %in% years)
levels(index.raw$gear) <- c('ATS2', 'ATS1', 'BTS')
aic1.1 <- Opt$AIC
out1.1 <- get.index.tmp('1: base case')
res1.1 <- get.resids.tmp('1: base case')
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
## Rerun with combined strata turned on
control$combinedoff <- TRUE; control$make_plots <- FALSE
savedir <- paste0(getwd(), '/test1_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic1.2 <- Opt$AIC
out1.2 <- get.index.tmp('1: base case')
res1.2 <- get.resids.tmp('1: base case')

## Add beta2 estimates
control$combinedoff <- FALSE; control$beta2temporal=TRUE
savedir <- paste0(getwd(), '/test2_combined')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic2.1 <- Opt$AIC
out2.1 <- get.index.tmp("2: + beta2 FE")
res2.1 <- get.resids.tmp("2: + beta2 FE")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test2_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic2.2 <- Opt$AIC
out2.2 <- get.index.tmp("2: + beta2 FE")
res2.2 <- get.resids.tmp("2: + beta2 FE")

## Add spatial component
control$combinedoff <- FALSE; control$n_omega1 <- control$n_omega2 <- 'IID'
savedir <- paste0(getwd(), '/test3_combined')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=0)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic3.1 <- Opt$AIC
out3.1 <- get.index.tmp("3: + spatial IID")
res3.1 <- get.resids.tmp("3: + spatial IID")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test3_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic3.2 <- Opt$AIC
out3.2 <- get.index.tmp("3: + spatial IID")
res3.2 <- get.resids.tmp("3: + spatial IID")


## Add spatiotemporal
control$combinedoff <- FALSE; control$n_eps1 <- control$n_eps2 <- 'IID'
savedir <- paste0(getwd(), '/test4_combined')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic4.1 <- Opt$AIC
out4.1 <- get.index.tmp("4: + spatiotemporal IID")
res4.1 <- get.resids.tmp("4: + spatiotemporal IID")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test4_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
results <- process.results(Opt, Obj, Inputs, model, space, savedir)
plot.vastfit(results)
aic4.2 <- Opt$AIC
out4.2 <- get.index.tmp("4: + spatiotemporal IID")
res4.2 <- get.resids.tmp("4: + spatiotemporal IID")



indices.all <- rbind(out1.1, out1.2, out2.1, out2.2, out3.1, out3.2,
                     out4.1, out4.2)
g <- ggplot(indices.all, aes(year, logdensity, color=type, linetype=model)) +
  geom_line(lwd=1, alpha=.8) + facet_grid(case~gear) + theme_bw()
ggsave('plots/index_buildup.png', g, width=8, height=9)
res.all <- rbind(res1.1, res1.2, res2.1, res2.2, res3.1, res3.2,
                     res4.1, res4.2)
g <- ggplot(subset(res.all, year==1), aes(obs, predicted, color=type)) + geom_point(alpha=.5) +
  facet_grid(case~gear) + geom_abline(slope=1, intercept=0)
ggsave('plots/index_buildup_resids.png', g, width=8, height=9)

aic.table <- cbind(c(aic1.1, aic2.1, aic3.1, aic4.1, aic5.1, aic6.1),
                   c(aic1.2, aic2.2, aic3.2, aic4.2, aic5.2, aic6.2))
delta.aic.table <- aic.table-apply(aic.table, 2, min)
delta.aic.table
