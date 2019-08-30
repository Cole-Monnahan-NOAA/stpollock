## A series of sensitivity analyses to run

## Different assumptions about the decorrelation range
source("run_kappa.R")
## Aniso turned off (isotropic) or fixed at the values from pcod
source("run_aniso.R")
## Test different spatial configurations: no space (NS), space only (S), or
## full spatiotemporal (ST)
source("run_spatialconfig.R")

source("run_temporal.R")


