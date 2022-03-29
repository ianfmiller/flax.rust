library(mgcv)
all_transects<-readRDS("~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_transects.RDS")
all_populations<-readRDS("~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_populations.RDS")
all_transects$landcover<-as.factor(all_transects$landcover)
all_transects$transect<-as.factor(all_transects$transect)

# landcover codes
##1=evergreen trees and shrubs
##2=deciduous trees greater than 2m tall
##3=meadow, grassland and subshrub
##4=persistent open water
##5=persistent snow and ice
##6=rock, bare soil, and sparse vegetation
##7=building or structure
##8=paved or other impervious surface
##9=irrigated pasture and other cultivated lands
##10=deciduous shrubs up to 2m tall
##11=evergreen forest understory and small gap
##12=deciduous forest understory and small gap

# where is flax?
mod0<-gam(flax.presence~
            s(elevation)+
            s(slope)+
            s(slope_southness)+
            s(slope_westness)+
            landcover,
          select=T,link=binomial(family="logit"),data=all_transects)
plot(mod0,pages=1)
summary(mod0)
## flax distribuion peaks around 3k meters
## flax likes steep slopes, to a point
## flax likes south facing slopes
## flax doesn't care about westness so long as it doesn't point it too far away from south
## flax likes (relative to evergreen trees and shrubs, in order of magnitude): 
### (3) meadow/grassland/shrubrub, 
### (6) rock/bare soil w/ sparse vegitation, irrigated pasture, 
### (9) irrigated pasture/cultivated
### (4) persistant open water???,
### (10) deciduous shrubs

## flax doesn't like (evergreen trees and shrubs, in order of magnitude)
### (8) paved or other impervious surface
### (2) deciduous trees greater than 2m tall

# where are flax populations densest?
mod1<-gam(density~
            s(elevation)+
            s(slope)+
            s(slope.southness)+
            s(slope.westness)+
            factor(mode.landcover),
          select=T,link=binomial(family="logit"),data=all_populations)
plot(mod1,pages=1)
summary(mod1)
### no relationship to anything

# where is flax rust, relative to where flax is
all_dis_transects<-all_transects[which(all_transects$flax.presence==1),]

mod2<-gam(incidence~
            s(elevation)+
            s(slope)+
            s(slope_southness)+
            s(slope_westness)+
            landcover,
          select=T,link=binomial(family="logit"),data=all_dis_transects)
plot(mod2,pages=1)
summary(mod2)

## no clear patterns w/ elevation, southness, westness
## steep slopes seem to have more disease
## flax rust likes (order of decreasing magnitude)
### (3) meadow/grassland/shrubrub
### (10) deciduous shrubs
### (6) rock/bare soil w/ sparse vegitation, irrigated pasture

#where is flax rust incidence highest amongst all pops?

## interpreting incidence as whether there was disease in the sub population
mod3a<-gam(ifelse(all_populations$num.D>0,1,0)~
            s(density)+
            s(elevation)+
            s(slope)+
            s(slope.southness)+
            s(slope.westness)+
            factor(mode.landcover),
          select=T,link=binomial(family="logit"),data=all_populations)
plot(mod3a,pages=1)
summary(mod3a)

### incidence higher for denser populations
### incidence higher for steeper slopes

## interpreting incidence as whether there was disease in the 'chunk'

mod3b<-gam(incidence~
             s(density)+
             s(elevation)+
             s(slope)+
             s(slope.southness)+
             s(slope.westness)+
             factor(mode.landcover),
           select=T,link=binomial(family="logit"),data=all_populations)
plot(mod3b,pages=1)
summary(mod3b)

### incidence higher for denser populations
### incidence higher for steeper slopes

#where is flax rust prevalence highest amongst all populations?
mod5<-gam(all_populations$num.D/(all_populations$num.H+all_populations$num.D)~
            s(density)+
            s(elevation)+
            s(slope)+
            s(slope.southness)+
            s(slope.westness)+
            factor(mode.landcover),
          select=T,link=binomial(family="logit"),data=all_populations)
plot(mod5,pages=1)
summary(mod5)

## nothing significant

#where is flax rust prevalence highest amongst diseased pops?
## had to drop some slope predictors due to insufficient data size
all_dis_populations<-all_populations[which(all_populations$incidence==1),]
all_dis_populations$prevalence<-all_dis_populations$num.D/(all_dis_populations$num.D+all_dis_populations$num.H)
mod6<-gam(prevalence~
            s(density)+
            s(elevation)+
            s(slope)+
            #s(slope.southness),
            #s(slope.westness),
            factor(mode.landcover),
          select=T,link=binomial(family="logit"),data=all_dis_populations)
plot(mod6,pages=1)
summary(mod6)

## highest for flat or steep slopes
## highest in meadows
