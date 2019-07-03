source("startup.R")
model <- 'combined'

### Test effect of increasing number of knots. Note I commented out the
### prepare.inputs line where I run teh model once for this test
nlls0 <- nlls <- inits0 <- inits <- results0 <- fits <- list()
nn <- seq(10, 50, by=5)
for(n in nn){
control <- list(beta2temporal=TRUE, n_x=n,
                ## n_eps1='IID', n_eps2='IID',
                n_omega2='IID', n_omega1="IID",
                beta1temporal=TRUE, filteryears=TRUE, combinedoff=TRUE,
                kappaoff=12, temporal=0, fixlambda=12, make_plots=FALSE)
savedir <- paste0(getwd(), '/test_', control$n_x)
source("prepare_inputs.R")
inits0[[n]] <- Obj$report(Obj$env$last.par)$jnll
nlls0[[n]] <- Obj$report()$jnll_comp
inits[[n]] <- Obj$fn()
Opt <- Optimize(obj=Obj, lower=TmbList$Lower, loopnum=3, getsd=FALSE,
                upper=TmbList$Upper,   savedir=savedir,
                newtonsteps=1, control=list(trace=10))
## results0[[n]] <- process.results(Opt, Obj, Inputs, model, space, savedir)
## plot.vastfit(results[[n]])
fits[[n]] <- Opt
nlls[[n]] <- Obj$report()$jnll_comp
}

par(mfrow=c(1,3))
plot(nn, unlist(inits0), type='b')
do.call(rbind, nlls0)
rowSums(do.call(rbind, nlls0))-unlist(inits0)
plot(nn, unlist(inits), type='b')
do.call(rbind, nlls)
rowSums(do.call(rbind, nlls))-unlist(inits)
plot(nn, unlist(sapply(fits, function(x) x$objective)), type='b',
     xlab='n_knots', y='NLL')
unlist(sapply(fits, function(x) x$max_grad))

par(mfrow=c(1,1))
tmp <- unlist(sapply(fits, function(x) x$objective))
plot(x=nn, type='b', xlab='n_knots', ylab='NLL',
     y=tmp-min(tmp))


## Look at why betas and lambdas can go so crazy
chains <- 1
options(mc.cores = \chains)
source('startup.R')
model <- 'combined'
control <- list(seed=121, beta2temporal=FALSE, filterdata=TRUE,
                n_eps1=0, n_eps2=0, n_omega1=0, n_omega2=0,
                beta1temporal=TRUE, combinedoff=FALSE)
savedir <- paste0(getwd(), '/mcmc_test')
source("prepare_inputs.R")
fit <- tmbstan(Obj, lower=TmbList$Lower, upper=TmbList$Upper, chains=chains,
               iter=500, open_progress=FALSE, warmup=200,
               init='last.par.best',
               control=list(max_treedepth=7))

df <- as.data.frame(fit)
plot(df[,1])
pairs(df[c(1,4,5,8,11)], gap=0, upper.panel=NULL)

jnll_comp <- list()
for(i in 1:nrow(df)){
  tmp <- Obj$report(df[i,-ncol(df)])
  jnll_comp[[i]] <- tmp$jnll_comp
}

jnll <- do.call(rbind, jnll_comp)[,11:12]
## jnllc <- apply(jnll, 2, function(x) x- mean(x))
## matplot(jnllc[1:100,], type='l')
plot(jnll)

par1 <- df[which.min(df[,1]),-11]
par2 <- df[which.max(df[,1]),-11]
data.frame(t(rbind(par1,par2)))
rep1 <- Obj$report(par1)
rep2 <- Obj$report(par2)


par(mfcol=c(2,3), mgp=c(1.5,.5,0), mar=c(3,3,.5,.5))
p1low <- rowSums(rep1$P1_iz)
p2low <- rowSums(rep1$P2_iz)
p1high <- rowSums(rep2$P1_iz)
p2high <- rowSums(rep2$P2_iz)
plot(p1low, ylab='P1', xlab='Data Row', ylim=range(c(p1low, p1high)))
points(p1high, col=2)
plot(p2low, ylab='P2', xlab='Data Row', ylim=range(c(p2low, p2high)))
points(p2high, col=2)
plot(rep1$R1_i, ylab='R1', xlab='Data Row', ylim=range(c(rep1$R1_i, rep2$R1_i)))
points(rep2$R1_i, col=2)
plot(rep1$R2_i, ylab='R2', xlab='Data Row', ylim=range(c(rep1$R2_i, rep2$R2_i)))
points(rep2$R2_i, col=2)
plot(rep1$LogProb1_i, ylab='LogProb1', xlab='Data Row',
     ylim=range(c(rep1$LogProb1_i, rep2$LogProb1_i)))
points(rep2$LogProb1_i, col=2)
legend("bottomright", legend=c('low', 'high'), col=c(1,2), pch=1)
plot(rep1$LogProb2_i, ylab='LogProb2', xlab='Data Row',
     ylim=range(c(rep1$LogProb2_i, rep2$LogProb2_)))
points(rep2$LogProb2_i, col=2)

plot(rep1$R1_i*rep1$R2_i, ylab='Density (R1*R2)', xlab='Data Row',
     ylim=range(c(rep1$R1_i*rep1$R2_i, rep2$R1_i*rep2$R2_i)))
points(rep2$R1_i*rep2$R2_i, col=2)


## Manually recreate calcs in R as a check. Do for one observation of each
## gear types

rep1$jnll_beta1
rep1$jnll_beta2
rep1$beta1_ft
rep1$beta1_tf
rep1$beta1_tc
beta1 <- rep1$beta1_tc[1,]
lambda1 <- par1$lambda1_k
lambda2 <- par1$lambda2_k

ind <- c(1,2000,6000)
rep1$P1_iz[ind,]
p1 <- matrix(-Inf, nrow=3, ncol=2)
p1[1,1] <- beta1[1] + lambda1
p1[1,2] <- beta1[2]+ lambda1
p1[2,1] <-  beta1[2]
p1[3,1] <-  beta1[3]
r1 <- 1-exp(-rowSums(exp(p1)))
r1
unique(rep1$R1_i)

R1 <- 1-exp(-exp(p1low))
plot(rep1$R1_i, R1, xlim=c(0,1), ylim=c(0,1)) ;abline(0,1)
plot(p1low, rep1$R1_i)
plot(p1low, R1)
