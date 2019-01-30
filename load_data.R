bts <- read.csv('data/bts.csv')
ats <- read.csv('data/ats.csv')
## The ats data is really high resolution so truncating this for now to
## make things faster and fix the mesh which is overly weighted to the ats
## data otherwise
ats <- ats[seq(1, nrow(ats), len=nrow(bts)),]
## fake depth data to test
bts$depth <- rnorm(nrow(bts))
ats$depth <- rnorm(nrow(ats))
