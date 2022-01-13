# Loading libraries: 
library(mgcv)
library(dplyr)
library(ggplot2)
library(ggpubr)


# Load data for each sampled point in "ribbon" along transects
all_transects<-readRDS("~/Desktop/Spatial Datasets/all_transects.RDS")
all_transects <- all_transects %>% mutate(landcover=factor(landcover))
str(all_transects) 

# Load data for each flax population
all_populations<-readRDS("~/Desktop/Spatial Datasets/all_populations.RDS")
all_populations<-subset(all_populations,density>0)
all_populations <- all_populations %>% mutate(prevalence = num.D/density) 



### Question: How does Elevation Affect Flax Locations and Disease Incidence/Prevalence 

    'Preliminary answer: elevation seems to affect flax locations but not necessarily disease dynamics'

    # Flax is differentially distributed at different elevations: 
    ggplot(data=all_transects, mapping = aes(x=elevation)) + 
        geom_point(aes(y = flax.presence)) + 
        labs(x= 'Elevation', y = 'Flax Presence', title = 'Elevation Effect on Flax Presence') 
     
    mod0<-gam(flax.presence~s(elevation),link=binomial,data=all_transects)
    summary(mod0)
    plot(mod0)
    
    #No linear effects of elevation on population density (kind of an unexpected result actually) 
    ggplot(data=all_populations, mapping = aes(x=elevation, y=density)) + 
        geom_point() + 
        geom_smooth(method='loess', se = F, color = 'red') + 
        stat_regline_equation(label.y = 125, label.x=3500,  aes(label = ..eq.label..)) +
        stat_regline_equation(label.y = 115, label.x=3500, aes(label = ..rr.label..))

    mod01 <- gam(density~s(elevation),data=all_populations)
    summary(mod01)
    plot(mod01)
    
    # Marginally different disease incidence at high elevations
    ggplot(data=all_populations, mapping = aes(x=elevation)) + 
        geom_point(aes(color = factor(incidence), y = num.H + num.D)) + 
        labs(x= 'Elevation', y = 'Number of Flax Plants', title = "Elevation's Effect on Incidence") + 
        geom_smooth(method='loess', se = F, aes(y = num.H + num.D, color = factor(incidence))) + 
        geom_smooth(method='loess', se = F, aes(y = num.H + num.D), color = 'black')

    # Upon further review, any effects of elevation on incidence are likely relatively small:
    incidence_mod <- glm(data=all_populations, incidence ~ elevation)
    summary(incidence_mod)
    incidence_mod1 <- gam(data=all_populations, incidence ~ s(elevation))
    plot(incidence_mod1)
    
    segmented_elevation <- cut(all_populations$elevation, seq(from=2600, to=3800, by=100))
    df <- as.data.frame(table(segmented_elevation))
    ggplot(data=all_populations, mapping=aes(x=factor(segmented_elevation), fill=factor(incidence))) + 
        labs(y='proportion', x='elevation groups', title = 'Ratios of disease incidence relative to elevation') + 
        geom_bar(position='fill') + 
        geom_text(data = df, aes(x=segmented_elevation, y=0.05, label=Freq), 
                  size=1.5, colour="black", inherit.aes=FALSE)
    
        # Check for the method of cutting then filling by a factor: 
        all_populations %>% filter(elevation<=2800 & elevation >= 2700) %>% pull(incidence) %>% length()
        all_populations %>% filter(elevation<=2800 & elevation >= 2700) %>% select(incidence) %>% table()
    
    
    # No significant effect of elevation on prevalence 
    prevalence_mod <- lm(data=all_populations, (num.D)/(num.H + num.D) ~ elevation) 
    summary(prevalence_mod)
    plot(all_populations$elevation, (all_populations$num.D)/(all_populations$num.H + all_populations$num.D))

    prevalence_mod1 <- gam((num.D/density)~s(elevation),data=all_populations)
    summary(prevalence_mod1)
    plot(prevalence_mod1)

    # Effects of elevation on metapopulation structure: 
    
    ggplot(data=all_populations, mapping=aes(x=elevation, y=nearest.pop.dist)) + 
        geom_point() + 
        geom_smooth(method='loess', se = F) + 
        labs(x= 'Elevation', y = 'Nearest Flax Population', title = "Elevation's Effect on Metapopulation Structure") 
    
        # Outlier seems to be influencing right tail of the graph 
   
    ggplot(data=all_populations, mapping=aes(x=elevation, y=nearest.D.pop.dist)) + 
        geom_point() + 
        geom_smooth(method='loess', se = F) + 
        labs(x= 'Elevation', y = 'Nearest Flax Population', title = "Elevation's Effect on Metapopulation Disease Structure") 
    
    # Significant result for metapopulation structures!
    mp_mod <- gam(nearest.pop.dist~s(elevation),data=all_populations)
    summary(mp_mod)
    plot(mp_mod)
    
    mp_mod2 <- gam(nearest.D.pop.dist~s(elevation),data=all_populations)
    summary(mp_mod2)
    plot(mp_mod2)
    
    
### Question: How does size of a diseased population affect local transmission?
    
    # Inconclusive, but seems like there is nothing interesting going on here 
    
    # flawed method b/c would double count but good for a first pass
    local_transmission <- all_populations %>% filter(incidence == 1 & nearest.pop.dist == nearest.D.pop.dist) 
    avg_size_transmitted <- mean(local_transmission$density)
    avg_size_D_transmitted <- mean(local_transmission$num.D)
    
    all_disease <- all_populations %>% filter(incidence == 1)
    avg_size <- mean(all_disease$density)
    avg_size_D <- mean(all_disease$num.D)
    
    # Total population size = insignificant 
    t.test(x = local_transmission$density, y = all_disease$density, 
           alternative= 'two.sided', 
           conf.level = 0.95)
    
    # Number of diseased = insignificant  
    t.test(x = local_transmission$num.D, y = all_disease$num.D, 
           alternative= 'two.sided', 
           conf.level = 0.95)
    
    # Density doesn't affect prevalence
    mod000<-lm(num.D/(num.D+num.H)~density,data=all_populations)
    summary(mod000)
           

    
### Question: How does connectivity relate to flax populations and disease?
    
    #Summary: Disease spreads from a population to its nearest neighbors
    
    #Distribution of connectivities: 
    factored_pop.dist <- cut(all_populations$nearest.pop.dist, seq(from=0, to=6000, by = 100))
    ggplot(data=all_populations, mapping = aes(x= factored_pop.dist)) + 
        geom_bar()
    
    #No real relationship between population density and connectivity: 
    ggplot(data=all_populations, mapping = aes(x= density, y=nearest.pop.dist)) + 
        geom_point() + 
        geom_smooth(method='loess', se = T) + 
        stat_regline_equation(label.y = 5250, label.x=100,  aes(label = ..eq.label..)) +
        stat_regline_equation(label.y = 5000, label.x=100, aes(label = ..rr.label..))
    
    
    # It seems that disease does seem to spread from one population to its nearest neighbor
    disease_nn <- (all_populations %>% filter(incidence == 1) %>% select(nearest.pop.dist)) == 
        (all_populations %>% filter(incidence == 1) %>% select(nearest.D.pop.dist))
    all_nn <- all_populations %>% select(nearest.pop.dist) == 
        all_populations %>% select(nearest.D.pop.dist)
    t.test(x = disease_nn, y = all_nn, 
           alternative= 'two.sided', 
           conf.level = 0.95)    
    
    # Significant at a 10% interval: 
    dc_mod <-glm(incidence~nearest.pop.dist,data=all_populations)
    summary(dc_mod)
    
   # no large relationship between connectivity and landcover 
    factored_pop.dist <- cut(all_populations$nearest.pop.dist, seq(from=0, to=6000, by = 600))
    ggplot(data= all_populations, mapping=aes(x=factor(mode.landcover), fill = factor(factored_pop.dist))) + 
        geom_bar(position='fill') + 
        labs(title = 'Landcover effect on Connectivity', x = 'Landcover', y = 'Proportion')
   
    cl_mod <- lm(nearest.pop.dist ~ factor(mode.landcover), data=all_populations)
    summary(cl_mod)
    
    
### Question: How does landcover affect flax and flax disease dynamics:  
    
    'Preliminary Conclusion: there are definitely landcover regions that contain more flax, but
    disease dynamics do not appear to be super related to landcover'
    
        'given that landcover 3 seems to be the most prevalent, maybe we could research whether climate change is expected to 
        increase or decrease the range of this landcover?'
        
    
    # It seems that flax is definitely found more in certain of the landcover categories: 
    df1 <- as.data.frame(table(all_transects$landcover))
    ggplot(data= all_transects, mapping=aes(x=factor(landcover), fill = factor(flax.presence))) + 
        geom_bar(position='fill') + 
        labs(title = 'Distribution of flax based on landcover', x = 'Landcover', y = 'Proprtion') + 
        geom_text(data = df1, aes(x=Var1, y=0.05, label=Freq), 
                  size=1.5, colour="black", inherit.aes=FALSE)
    
    ggplot(data= all_populations, mapping=aes(x=factor(mode.landcover)))+ 
        geom_bar() 
    
    
    # Landcover x density 
    factored_pop.den <- cut(all_populations$density, seq(from=0, to=130, by = 13))
    ggplot(data= all_populations, mapping=aes(x=factor(factored_pop.den), fill = factor(mode.landcover))) + 
        geom_bar() + 
        labs(title = 'Density x Landcover', x = 'Density')
    
    
    # Flax rust incidence just seems to be higher in the landcover categories where there is flax (duh): 
    ggplot(data= all_populations, mapping=aes(x=factor(mode.landcover), fill = factor(incidence))) + 
        geom_bar(position='fill') + 
        labs(title = 'Landcover effect on incidence', x = 'Landcover', y = 'Proportion')
    
    mod7<-glm(data=all_populations, incidence~as.factor(mode.landcover))
    summary(mod7)
    
    # Same as with incidence: prevalence is just higher where there is more flax 
    prevalence_factor <- factor(cut(all_populations$prevalence, seq(from = 0, to=1, by=0.1), right = F)) 
    ggplot(data= all_populations, mapping=aes(x=factor(mode.landcover), fill = prevalence_factor)) + 
        geom_bar(position='fill') + 
        labs(title = 'Landcover effect on prevalence', x = 'Landcover', y= "Proprtion of prevalence category")
    
    mod8<-lm(num.D/(num.D+num.H)~as.factor(mode.landcover)+elevation,data=all_populations)
    summary(mod8)
    
    
# Bonus Transect Analysis: 
    
    #Landcover: 
    ggplot(data=all_populations, mapping=aes(x=factor(transect), fill = factor(mode.landcover))) + 
        geom_bar(position='fill') 
    
    #Density: 
    tr_dn<- aggregate(all_populations$density, by=list(all_populations$transect), FUN=mean) %>% rename(values = x) %>% group_by(values) %>% arrange(.by_group =T)

    #Elevation:
    tr_ev<- aggregate(all_populations$elevation, by=list(all_populations$transect), FUN=mean) %>% rename(values = x) %>% group_by(values) %>% arrange(.by_group =T)
    
  
  