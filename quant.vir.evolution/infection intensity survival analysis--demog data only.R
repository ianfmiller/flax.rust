setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")

demog<-subset(demog,status %in% c("H","D","X"))
demog<-subset(demog,final.status %in% c("H","D","X"))

tag<-c()
status<-c()
death<-c()
height<-c()
ninfl<-c()
nDinfl<-c()
pDinfl<-c()

for(index in unique(demog$tag))
{
  
  stat1<-NA
  stat2<-NA
  death1<-NA
  death2<-NA
  height1<-NA
  height2<-NA
  ninfl1<-NA
  ninfl2<-NA
  nDinfl1<-NA
  nDinfl2<-NA
  
  sub.data<-subset(demog,tag==index)
  if(all("2018" %in% sub.data$year,"2019" %in% sub.data$year))
  {
    stat1<-as.character(sub.data[which(sub.data$year=="2018"),"final.status"])
    death1<-ifelse(sub.data[which(sub.data$year=="2019"),"status"]=="X",1,0)
    height1<-as.numeric(sub.data[which(sub.data$year=="2018"),"height.cm"])
    ninfl1<-as.numeric(sub.data[which(sub.data$year=="2018"),"num.infl"])
    nDinfl1<-as.numeric(sub.data[which(sub.data$year=="2018"),"num.D.infl"])
  }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
  {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    death2<-ifelse(sub.data[which(sub.data$year=="2020"),"status"]=="X",1,0)
    height2<-as.numeric(sub.data[which(sub.data$year=="2019"),"height.cm"])
    ninfl2<-as.numeric(sub.data[which(sub.data$year=="2019"),"num.infl"])
    nDinfl2<-as.numeric(sub.data[which(sub.data$year=="2019"),"num.D.infl"])
  }
  
  if(!any(is.na(c(stat1,death1)))) {tag<-c(tag,index);status<-c(status,stat1);death<-c(death,death1);height<-c(height,height1);ninfl<-c(ninfl,ninfl1);nDinfl<-c(nDinfl,nDinfl1)}
  if(!any(is.na(c(stat2,death2)))) {tag<-c(tag,index);status<-c(status,stat2);death<-c(death,death2);height<-c(height,height2);ninfl<-c(ninfl,ninfl2);nDinfl<-c(nDinfl,nDinfl2)}
  
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,status=status,death=death,height=height,ninfl=ninfl,nDinfl=nDinfl,pDinfl=nDinfl/ninfl)
surv.data<-droplevels(surv.data)
surv.data2<-surv.data[-which(is.na(surv.data$pDinfl)),]

#mod1<-glm(death~pDinfl*height*ninfl,data=surv.data2,family=binomial(link="logit"))
#mod2<-glm(death~pDinfl+height*ninfl,data=surv.data2,family=binomial(link="logit"))
#mod3<-glm(death~pDinfl*height+ninfl,data=surv.data2,family=binomial(link="logit"))
#mod4<-glm(death~pDinfl*ninfl+height,data=surv.data2,family=binomial(link="logit"))
#mod5<-glm(death~pDinfl+height+ninfl,data=surv.data2,family=binomial(link="logit"))
#mod6<-glm(death~pDinfl*height,data=surv.data2,family=binomial(link="logit"))
#mod7<-glm(death~pDinfl+height,data=surv.data2,family=binomial(link="logit"))
#mod8<-glm(death~pDinfl*ninfl,data=surv.data2,family=binomial(link="logit"))
#mod9<-glm(death~pDinfl+ninfl,data=surv.data2,family=binomial(link="logit"))
#mod10<-glm(death~height*ninfl,data=surv.data2,family=binomial(link="logit"))
#mod11<-glm(death~height+ninfl,data=surv.data2,family=binomial(link="logit"))
#mod12<-glm(death~pDinfl,data=surv.data2,family=binomial(link="logit"))
#mod13<-glm(death~height,data=surv.data2,family=binomial(link="logit"))
#mod14<-glm(death~ninfl,data=surv.data2,family=binomial(link="logit"))
#AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)
final.mod<-glm(death~pDinfl+height,data=surv.data2,family=binomial(link="logit"))

ilogit<-function(x) {exp(x)/(exp(x)+1)} 

plot(surv.data2$pDinfl,surv.data2$death,xlab="% stems infected",ylab="prob.death")
curve(ilogit(coef(final.mod)[1]+coef(final.mod)[2]*x+coef(final.mod)[3]*quantile(surv.data2$height,.05)),add=T,col=heat.colors(5)[1])
curve(ilogit(coef(final.mod)[1]+coef(final.mod)[2]*x+coef(final.mod)[3]*quantile(surv.data2$height,.25)),add=T,col=heat.colors(5)[2])
curve(ilogit(coef(final.mod)[1]+coef(final.mod)[2]*x+coef(final.mod)[3]*quantile(surv.data2$height,.5)),add=T,col=heat.colors(5)[3])
curve(ilogit(coef(final.mod)[1]+coef(final.mod)[2]*x+coef(final.mod)[3]*quantile(surv.data2$height,.75)),add=T,col=heat.colors(5)[4])
curve(ilogit(coef(final.mod)[1]+coef(final.mod)[2]*x+coef(final.mod)[3]*quantile(surv.data2$height,.95)),add=T,col=heat.colors(5)[5])
legend("topright",legend=c("5%","25%","50%","75%","95%"),lwd=1,col=heat.colors(5),title = "quantile height")


#mod1<-glm(death~nDinfl*height*ninfl,data=surv.data,family=binomial(link="logit"))
#mod2<-glm(death~nDinfl+height*ninfl,data=surv.data,family=binomial(link="logit"))
#mod3<-glm(death~nDinfl*height+ninfl,data=surv.data,family=binomial(link="logit"))
#mod4<-glm(death~nDinfl*ninfl+height,data=surv.data,family=binomial(link="logit"))
#mod5<-glm(death~nDinfl+height+ninfl,data=surv.data,family=binomial(link="logit"))
#mod6<-glm(death~nDinfl*height,data=surv.data,family=binomial(link="logit"))
#mod7<-glm(death~nDinfl+height,data=surv.data,family=binomial(link="logit"))
#mod8<-glm(death~nDinfl*ninfl,data=surv.data,family=binomial(link="logit"))
#mod9<-glm(death~nDinfl+ninfl,data=surv.data,family=binomial(link="logit"))
#mod10<-glm(death~height*ninfl,data=surv.data,family=binomial(link="logit"))
#mod11<-glm(death~height+ninfl,data=surv.data,family=binomial(link="logit"))
#mod12<-glm(death~nDinfl,data=surv.data,family=binomial(link="logit"))
#mod13<-glm(death~height,data=surv.data,family=binomial(link="logit"))
#mod14<-glm(death~ninfl,data=surv.data,family=binomial(link="logit"))
#AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)

final.mod2<-glm(death~nDinfl+ninfl,data=surv.data,family=binomial(link="logit"))

ilogit<-function(x) {exp(x)/(exp(x)+1)} 

plot(surv.data$nDinfl,surv.data$death,xlab="N stems infected",ylab="prob.death")
curve(ilogit(coef(final.mod2)[1]+coef(final.mod2)[2]*x+coef(final.mod2)[3]*quantile(surv.data$ninfl,.05,na.rm = T)),add=T,col=heat.colors(5)[1])
curve(ilogit(coef(final.mod2)[1]+coef(final.mod2)[2]*x+coef(final.mod2)[3]*quantile(surv.data$ninfl,.25,na.rm = T)),add=T,col=heat.colors(5)[2])
curve(ilogit(coef(final.mod2)[1]+coef(final.mod2)[2]*x+coef(final.mod2)[3]*quantile(surv.data$ninfl,.5,na.rm = T)),add=T,col=heat.colors(5)[3])
curve(ilogit(coef(final.mod)[1]+coef(final.mod2)[2]*x+coef(final.mod2)[3]*quantile(surv.data$ninfl,.75,na.rm = T)),add=T,col=heat.colors(5)[4])
curve(ilogit(coef(final.mod2)[1]+coef(final.mod2)[2]*x+coef(final.mod2)[3]*quantile(surv.data$ninfl,.95,na.rm = T)),add=T,col=heat.colors(5)[5])
legend("topright",legend=c("5%","25%","50%","75%","95%"),lwd=1,col=heat.colors(5),title = "quantile N stems")

#AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)