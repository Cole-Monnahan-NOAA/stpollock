## Base case simplified model
control <- list(seed=121, beta2temporal=FALSE, beta1temporal=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                combinedoff=FALSE, filteryears=TRUE, make_plots=TRUE,
                kappaoff=0, temporal=0, fixlambda=12)
savedir <- paste0(getwd(), '/test1_combined')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
index.raw <- read.csv(paste0(savedir, '/index.raw.csv'))
index.raw <- subset(index.raw, type=='Naive Spatial')
levels(index.raw$gear) <- c('ATS2', 'ATS1', 'BTS')
out1.1 <- get.index.tmp('1: base case')
## Rerun with combined strata turned on
control$combinedoff <- TRUE; control$make_plots <- FALSE
savedir <- paste0(getwd(), '/test1_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out1.2 <- get.index.tmp('1: base case')

## Add beta2 estimates
control$combinedoff <- FALSE; control$beta2temporal=TRUE
savedir <- paste0(getwd(), '/test2_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out2.1 <- get.index.tmp("2: + beta2 FE")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test2_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out2.2 <- get.index.tmp("2: + beta2 FE")

## Add Omega1 IID
control$combinedoff <- FALSE; control$n_omega1='IID'
savedir <- paste0(getwd(), '/test3_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out3.1 <- get.index.tmp("3: + omega1 IID")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test3_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out3.2 <- get.index.tmp("3: + omega1 IID")

## Add Omega2 IID
control$combinedoff <- FALSE; control$n_omega2='IID'
savedir <- paste0(getwd(), '/test4_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out4.1 <- get.index.tmp("4: + omega2 IID")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test4_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out4.2 <- get.index.tmp("4: + omega2 IID")

## Add Omega1 and Omega2 L's
control$combinedoff <- FALSE; control$n_omega2=3; control$n_omega1=3
savedir <- paste0(getwd(), '/test5_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out5.1 <- get.index.tmp("5: + omegas SFA")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test5_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out5.2 <- get.index.tmp("5: + omegas SFA")

## Add IID eps1
control$combinedoff <- FALSE; control$n_eps1="IID"
savedir <- paste0(getwd(), '/test6_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out6.1 <- get.index.tmp("6: + eps1 IID")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test6_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out6.2 <- get.index.tmp("6: + eps1 IID")

## Add full rank epsilon1
control$combinedoff <- FALSE; control$n_eps1=3
savedir <- paste0(getwd(), '/test7_combine')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out7.1 <- get.index.tmp("7: + eps1 SFA")
control$combinedoff <- TRUE
savedir <- paste0(getwd(), '/test7_combinedoff')
source("prepare_inputs.R")
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=TRUE,
                upper=TmbList$Upper, savedir=savedir, newtonsteps=1)
out7.2 <- get.index.tmp("7: + eps1 SFA")

indices.all <- rbind(out1.1, out1.2, out2.1, out2.2, out3.1, out3.2,
                     out4.1, out4.2, out5.1, out5.2, out6.1, out6.2,
                     out7.1, out7.2)
g <- ggplot(indices.all, aes(year, logdensity, color=type, linetype=model)) +
  geom_line(lwd=1, alpha=.8) + facet_grid(case~gear) + theme_bw()
ggsave('plots/index_buildup.png', g, width=8, height=9)
