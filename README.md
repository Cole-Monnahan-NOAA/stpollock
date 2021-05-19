# stpollock 

A spatio-temporal CPUE analysis combining bottom and acoustic
trawl survey data for eastern Bering sea pollock


This repo contains code to reproduce the following paper:

Cole C Monnahan, James T Thorson, Stan Kotwicki, Nathan
Lauffenburger, James N Ianelli, Andre E Punt, Incorporating
vertical distribution in index standardization accounts for
spatiotemporal availability to acoustic and bottom trawl gear for
semi-pelagic species, ICES Journal of Marine Science, 2021;,
fsab085, https://doi.org/10.1093/icesjms/fsab085 

A demo of the method is posted at: 
https://github.com/James-Thorson-NOAA/VAST/wiki/Combine-acoustic-and-bottom-trawl-data

Note that this example uses maximum likelihood, while the paper
above uses Bayesian inference via tmbstan. See paper for a
discussion of why.

To reproduce the analysis, clone the repo, open the script
"run_analysis.R" and step through calculations. The whole thing
may take a month or so to run.

The result files are far too large to save so are
excluded. Otherwise all files necessary to run the analysis are
included. 

 
 
 
 
