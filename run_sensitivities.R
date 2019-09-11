## A series of sensitivity analyses to run. Each of these takes a
## long time to run, and specific settings can be found in each
## script. The result is a plot in the 'plots' folder and a
## result file in 'results' folder.

## Different assumptions about the decorrelation range
source("sensitivities/run_kappa.R")

## Aniso turned off (isotropic) or fixed at the values from
## previous run
source("sensitivities/run_aniso.R")

## Test different spatial configurations: no space (NS), space only (S), or
## full spatiotemporal (ST)
source("sensitivities/run_spatialconfig.R")

## source("sensitivities/run_temporal.R")

source("sensitivities/run_catchability.R")

