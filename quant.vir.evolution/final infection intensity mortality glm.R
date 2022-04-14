# load data
library(lme4)
setwd("~/Documents/GitHub/flax.rust/data")
intensity<-read.csv("Withinhost.csv")
demog<-read.csv("~/Downloads/Demography.csv") #git hub version doesn't have 2021 yet

# clean intensity data
intensity<-intensity[,c("Year","Site","Tag","Date","N.D.Stems","max.height","stem.height","percent.tissue.infected","length.tissue.infected","N.pustules.middle"),]
intensity<-intensity[which(intensity$Year %in% c("2019","2020")),]
intensity[which(intensity$Site=="BTDO"),"Site"]<-"BT"

## conditions for inclusion

### 2019
#### CC: 6/20 or 6/21 through 7/30
#### BT: 6/24 through 7/30
#### GM: not tracked whole season
#### HM: not tracked whole season
### 2020
#### CC: 6/25 through 7/28
#### BT 6/23 through 7/29
#### GM: 6/16 or 6/23 through 7/28
#### HM: 6/25 through 7/28

first.dates<-c(as.Date("2019/06/21"),as.Date("2019/06/24"),NA,NA,as.Date("2020/06/25"),as.Date("2020/06/23"),as.Date("2020/06/23"),as.Date("2020/06/25"))
last.dates<-c(as.Date("2019/07/30"),as.Date("2019/07/30"),NA,NA,as.Date("2020/07/27"),as.Date("2020/07/29"),as.Date("2020/07/28"),as.Date("2020/07/28"))

infection.intensities<-data.frame("year"=numeric(),"site"=character(),"tag"=character(),"infection.intensity"=numeric(),"height"=numeric())
i<-1
for(year in c("2019","2020"))
{
  for(site in c("CC","BT","GM","HM"))
  {
    sub.intensity<-intensity[which((intensity$Year==year) & (intensity$Site==site)),]
    sub.intensity.last.obs<-sub.intensity[which(as.Date(sub.intensity$Date,tryFormats = "%m/%d/%y")==last.dates[i]),]
    
    for(tag in unique(sub.intensity.last.obs$Tag))
    {
      sub.intensity.last.obs.tag<-sub.intensity.last.obs[which(sub.intensity.last.obs$Tag==tag),]
      if(!(sub.intensity.last.obs.tag[1,]$length.tissue.infected==""))
      {
        mean.lengths=mean(as.numeric(sub.intensity.last.obs.tag$length.tissue.infected),na.rm=T)
      } else
      {
        mean.lengths=mean(as.numeric(sub.intensity.last.obs.tag$stem.height)*as.numeric(sub.intensity.last.obs.tag$percent.tissue.infected),na.rm=T)
      }
      mean.n.pustules=mean(as.numeric(sub.intensity.last.obs.tag$N.pustules.middle),na.rm=T)
      N.D.Stems<-sub.intensity.last.obs.tag[1,"N.D.Stems"]
      inf.intens<-mean.lengths*mean.n.pustules*N.D.Stems
      height=as.numeric(sub.intensity.last.obs.tag[1,"max.height"])
      infection.intensities<-rbind(infection.intensities,data.frame("year"=year,"site"=site,"tag"=tag,"infection.intensity"=inf.intens,"height"=height))
    }
    i<-i+1
  }
}

infection.intensities<-infection.intensities[-which(infection.intensities$infection.intensity==0),]
infection.intensities<-infection.intensities[-which(is.na(infection.intensities$infection.intensity)),]

# join with mortality data

demog[which(demog$Site=="CCDO"),"Site"]<-"CC"
demog[which(demog$Site=="BTDO"),"Site"]<-"BT"
demog[which(demog$Site=="GMDO"),"Site"]<-"GM"

mortality<-c()
for(i in 1:nrow(infection.intensities))
{
  year<-as.numeric(infection.intensities[i,"year"])
  site<-infection.intensities[i,"site"]
  tag<-infection.intensities[i,"tag"]
  
  demog.sub<-demog[which(demog$year==(year+1)),]
  demog.sub<-demog.sub[which(demog.sub$Site==site),]
  demog.sub<-demog.sub[which(demog.sub$tag==tag),]
  
  mort<-NA
  if(nrow(demog.sub)>0)
  {
    if(demog.sub$status %in% c("D","H")) {mort<-0}
    if(demog.sub$status=="X") {mort<-1}
  }
  print(paste0("i = ",i," data = ",demog.sub$status," coded = ",mort))
  mortality<-c(mortality,mort)
}

data<-cbind(infection.intensities,"mortality"=mortality)

mod<-glmer(mortality~infection.intensity+height+(1|site)+(1|year),data=data,family = binomial("logit"))

