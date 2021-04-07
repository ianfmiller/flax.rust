#these scripts are for comparing the epi and within.host data sets to ID what data to include in the GLM analysis of inf outcome ~ FOI

# load + prep data
## load data
epi<-read.csv("~/Documents/GitHub/flax.rust/data/Epidemiology.csv")
within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
plant.locs<-read.csv("~/Documents/GitHub/flax.rust/data/plant locs.csv")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")

## prep data
epi<-epi[which(epi$Year==2020),]
epi$Date.First.Observed.Diseased<-as.Date(epi$Date.First.Observed.Diseased,tryFormats = "%m/%d/%Y")
within.host<-within.host[which(within.host$Year==2020),]
within.host$Date<-as.Date(within.host$Date,tryFormats = "%m/%d/%Y")

#visualize differences between epi and within host data sets
for (site in unique(epi$Site))
{
  sub.data.1a<-epi[which(epi$Site==site),]
  sub.data.1a<-sub.data.1a[order(sub.data.1a$Date.First.Observed.Diseased),]
  sub.data.1b<-within.host[which(within.host$Site==site),c("Year","Site","Tag","Date")]
  sub.data.1b<-unique(sub.data.1b)
  sub.data.1b<-sub.data.1b[order(sub.data.1b$Date),]
  for (date in unique(sub.data.1a$Date.First.Observed.Diseased))
  {
    sub.data.2a<-sub.data.1a[which(sub.data.1a$Date.First.Observed.Diseased<=date),]
    sub.data.2b<-sub.data.1b[which(sub.data.1b$Date==date),]
    if(dim(sub.data.2b)[1]>=1)
    {
      plot(0,0,xlim=c(0,10),ylim=c(0,20),type="n",main=paste0(site," date = ",as.Date(date,origin = "1970-01-01")))
      # data without observed plant inf intens
      sub.data.3a<-sub.data.2a[intersect(which(!(sub.data.2a$Tag %in% sub.data.2b$Tag)),which(!(sub.data.2a$notes=="seedling"))),]
      # seedling data to be assigned arbitrary starting plant inf intens value
      sub.data.3b<-sub.data.2a[intersect(which(!(sub.data.2a$Tag %in% sub.data.2b$Tag)),which(sub.data.2a$notes=="seedling")),]
      # data with observed plant inf intens
      sub.data.3c<-sub.data.2a[which(sub.data.2a$Tag %in% sub.data.2b$Tag),]
      points(sub.data.3a$X+sub.data.3a$x,sub.data.3a$Y+sub.data.3a$y,col="black") 
      points(sub.data.3b$X+sub.data.3b$x,sub.data.3b$Y+sub.data.3b$y,col="green") 
      points(sub.data.3c$X+sub.data.3c$x,sub.data.3c$Y+sub.data.3c$y,col="red") 
    }
    Sys.sleep(2)
  }
}

# compare data dimensions
for (site in unique(epi$Site))
{
  sub.data.1<-epi[which(epi$Site==site),]
  sub.data.1<-sub.data.1[order(sub.data.1$Date.First.Observed.Diseased),]
  plot(0,0,xlim=c(0,10),ylim=c(0,20),type="n",main=site)
  i<-0
  for (date in unique(sub.data.1$Date.First.Observed.Diseased))
  {
    sub.data.2<-sub.data.1[which(sub.data.1$Date.First.Observed.Diseased==date),]
    i<-i+dim(sub.data.2)[1]
    print(paste0("site = ",site," date = ",as.Date(date,origin = "1970-01-01")," n tot = ",i))
  }
}

for (site in unique(within.host$Site))
{
  sub.data.1<-within.host[which(within.host$Site==site),c("Year","Site","Tag","Date")]
  sub.data.1<-unique(sub.data.1)
  sub.data.1<-sub.data.1[order(sub.data.1$Date),]
  plot(0,0,xlim=c(0,10),ylim=c(0,20),type="n",main=site)
  for (date in unique(sub.data.1$Date))
  {
    sub.data.2<-sub.data.1[which(sub.data.1$Date==date),]
    print(paste0("site = ",site," date = ",date," n tot = ",dim(sub.data.2)[1]))
  }
}