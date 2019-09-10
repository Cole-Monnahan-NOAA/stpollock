## A series of sensitivity analyses to run

## Different assumptions about the decorrelation range
source("sensitivities/run_kappa.R")
## Aniso turned off (isotropic) or fixed at the values from pcod
source("sensitivities/run_aniso.R")
## Test different spatial configurations: no space (NS), space only (S), or
## full spatiotemporal (ST)
source("sensitivities/run_spatialconfig.R")

source("sensitivities/run_temporal.R")

source("sensitivities/run_catchability.R")

