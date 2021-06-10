# script to build dataset of observed/predicted plant heights on dates relevant for foi dataset building

## setup

### convenience params
data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09"))

### load data
corrected.epi<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.epi.RDS")
corrected.locs<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.locs.RDS")
within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
healthyplants<-read.csv("~/Documents/GitHub/flax.rust/data/healthyplants.csv")

## analysis 

### for plant in corrected.locs
### for date in data.dates at plant site
### if height data exists--store it
### if height data doesn't exist--hindcast/forecast