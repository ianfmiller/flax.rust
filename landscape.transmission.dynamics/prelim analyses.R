## density predicts incidence, not prevalence
flax_pops<-all_populations[which(all_populations$density>0),]
par(mfrow=c(1,2))
plot(flax_pops$density,flax_pops$incidence)
plot(flax_pops$density,flax_pops$num.D/(flax_pops$num.D+flax_pops$num.H))
mod1<-glm(incidence~density,data = flax_pops)
mod2<-lm(num.D/(num.D+num.H)~density,data=flax_pops)
summary(mod1)
summary(mod2)

## connectivity doesn't predict anything
plot(flax_pops$nearest.pop.dist,flax_pops$incidence)
plot(flax_pops$nearest.pop.dist,flax_pops$num.D/(flax_pops$num.D+flax_pops$num.H))
mod3<-glm(incidence~nearest.pop.dist+elevation,data=flax_pops)
mod4<-lm(num.D/(num.D+num.H)~density+elevation,data=flax_pops)
summary(mod3)
summary(mod4)

## landcover xxx flax presence

#plot(jitter(all_transects$landcover),jitter(all_transects$flax.presence))
mod5<-lm(flax.presence~as.factor(landcover)+elevation,data=all_transects)
summary(mod5)

## landcover xx density
plot(jitter(all_populations$mode.landcover),all_populations$density)
plot(all_populations$p.landcover.3,all_populations$density)
mod6<-lm(density~as.factor(mode.landcover)+elevation,data=all_populations)
mod7<-lm(density~p.landcover.3,data=all_populations)
summary(mod6)
summary(mod7)

## landcover incidence
plot(jitter(all_populations$mode.landcover),all_populations$incidence)
mod7<-glm(incidence~as.factor(mode.landcover)+elevation,data=all_populations)
summary(mod7)

## landcover prevalence
plot(jitter(all_populations$mode.landcover),all_populations$num.D/(all_populations$num.H+all_populations$num.D))
mod8<-lm(num.D/(num.D+num.H)~as.factor(mode.landcover)+elevation,data=all_populations)
summary(mod8)

