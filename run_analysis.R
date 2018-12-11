## Script to run the analysis for combining bottom trawl and acoustic data
## in a spatiotemporal index standardization.

## Started 8/2018 by Cole

### Step 1: Prepare workspace for analysis
source("startup.R")

### Step 2: Fit models to the real data
source("fit_models.R")

### Step 3: Run a simulation test
source("run_simulation.R")

### Step 3: Check model diagnostics
source("eval_models.R")
source("make_plots.R")
