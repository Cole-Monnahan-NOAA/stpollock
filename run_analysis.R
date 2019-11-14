## Script to run the analysis for combining bottom trawl and acoustic data
## in a spatiotemporal index standardization.

## Started 8/2018 by Cole

### Step 1: Prepare workspace for analysis
source("startup.R")

### Step 2: Fit models to the real data
source("fit_basecase.R")

### Step 3: Run a simulation test
source("run_simulation.R")

### Step 3: A series of sensitivity analyses to run. Each of
### these takes a long time to run, and specific settings can be
### found in each script. The result is a plot in the 'plots'
### folder and a result file in 'results' folder.
## Different assumptions about the decorrelation range
source("sensitivities/run_kappa.R")
## Aniso turned off (isotropic) or fixed at the values from
## previous run
source("sensitivities/run_aniso.R")
## Test different spatial configurations: no space (NS), space only (S), or
## full spatiotemporal (ST)
source("sensitivities/run_spatialconfig.R")
source("sensitivities/run_catchability.R")
source("sensitivities/run_resolution.R")
source("sensitivities/run_inflated0.R")

### Step 3: Check model diagnostics
source("eval_models.R")
source("make_plots.R")



