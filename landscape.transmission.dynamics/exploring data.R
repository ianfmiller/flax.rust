# exploring data
library(mgcv)

## load data for each sampled point in "ribbon" along transects
all_transects<-readRDS("~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_transects.RDS")
head(all_transects)

## load data for each flax population
all_populations<-readRDS("~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_populations.RDS")
all_populations<-subset(all_populations,density>0)
head(all_populations)

# some things to get you started...
plot(all_transects$elevation,all_transects$flax.presence)
mod0<-gam(flax.presence~s(elevation),link=binomial,data=all_transects)
plot(mod0)

plot(all_populations$density,all_populations$num.D/all_populations$density)

