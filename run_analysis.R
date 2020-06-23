## Script to run the analysis for combining bottom trawl and acoustic data
## in a spatiotemporal index standardization.

## Started 8/2018 by Cole

### Step 1: Prepare workspace for analysis
source("startup.R")
packageVersion('TMB')                   # 1.7.16
packageVersion('tmbstan')               # 1.0.2
packageVersion('rstan')                 # 2.19.3

### Step 2: Fit models to the pollock data
source("fit_basecase.R")

### Step 3: Run a simulation test
source("run_simulation.R")

### Step 4: A series of sensitivity analyses to run. Each of
### these takes a long time to run, and specific settings can be
### found in each script. The result is a plot in the 'plots'
### folder and a result file in 'results' folder. These all use
### n_x=200 (except spatial configuration). Note that the kappa=1
### case in run_kappa is used as a base case 200 run for the
### others so this needs to be run first.

## Different assumptions about the decorrelation range. Note
source("sensitivities/run_kappa.R")

## Aniso turned off (isotropic) or fixed at the values from
## previous run
source("sensitivities/run_aniso.R")

## Test different spatial configurations: no space (NS), space only (S), or
## full spatiotemporal (ST)
source("sensitivities/run_spatialconfig.R")

## Configurations of catchability
source("sensitivities/run_catchability.R")

## Spatial resolution (50,100,200,400)
source("sensitivities/run_resolution.R")

## Effect of inflated zeroes in SE corner which is unsampled by AT
source("sensitivities/run_inflated0.R")

## Effect of effective fishing height (3m vs 16m)
source("sensitivities/run_breakpoints.R")




