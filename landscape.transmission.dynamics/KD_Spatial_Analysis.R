### Flax Spatial Patterns of Disease:


## Environment Set Up: 
setwd("/Users/keenanduggal/desktop/Spatial Datasets")
library(dplyr)
library(ncf) 
library(sf)
library(sp)
library(ncf) 
library(spdep)

# Currently doesn't use Ian's Data Extraction Method, instead goes back to raw flax population data: 
flax_pops <- read.csv("landscape_transects.csv")
flax_pops$transect <- as.factor(flax_pops$transect)
flax_pops <- flax_pops %>% dplyr::select(transect, start.long, start.lat, num.H, num.D) %>% group_by(transect) %>% head(-1) 
flax_pops$subpop.incidence <- as.integer(flax_pops$num.D > 0)
flax_pops$subpop.prevalence <- (flax_pops$num.D / (flax_pops$num.D + flax_pops$num.H))



## Spatial Analysis Using All Populations / Transects Grouped Together: 


# Mantel Test: Tests whether populations spatially close to one another also have similar disease statuses:

Mantel_incidence = mantel.test(x = flax_pops$start.long, y = flax_pops$start.lat, z = flax_pops$subpop.incidence)
    # Very weak correlation that would  say further populations have more similar disease incidences

    #Different configuration of the mantel test (very similar results)
    n1 = length(flax_pops$subpop.incidence)
    mu1 = mean(flax_pops$subpop.incidence)
    sig1 = sd(flax_pops$subpop.incidence) * (n1 - 1)/n1
    zscale1 = (flax_pops$subpop.incidence - mu1)/sig1
    rho1 = outer(zscale1, zscale1)
    dst1 = as.matrix(dist(flax_pops[, c("start.long", "start.lat")]))
    mantel_incidence1 = mantel.test(M1 = rho1, M2 = dst1)
     # very weak correlation -- opposite of what we want 

Mantel_prevalence = mantel.test(x = flax_pops$start.long, y = flax_pops$start.lat, z = flax_pops$subpop.prevalence)
    # Weak correlation suggesting that further populations have less similar disease incidences 

    #Different configuration of the mantel test (very similar results) 
    n2 = length(flax_pops$subpop.prevalence)
    mu2 = mean(flax_pops$subpop.prevalence)
    sig2 = sd(flax_pops$subpop.prevalence) * (n2 - 1)/n2
    zscale2 = (flax_pops$subpop.prevalence - mu2)/sig2
    rho2 = outer(zscale2, zscale2)
    dst2 = as.matrix(dist(flax_pops[, c("start.long", "start.lat")]))
    mantel_prevalence2 = mantel.test(M1 = rho2, M2 = dst2)


# Correlogram: shows autocorrelation as a function of distance 
    correlog_prevtest <- correlog(x=flax_pops$start.long, y=flax_pops$start.lat, z=flax_pops$subpop.prevalence,increment = 0.025)
    plot(correlog_prevtest)
    correlog_inctest <- correlog(x=flax_pops$start.long, y=flax_pops$start.lat, z=flax_pops$subpop.incidence,increment = 0.025)
    plot(correlog_inctest)
    # no filled circles ==> no signficicance 

    
# Non parametric spatial correlation: 
    spline_prev =spline.correlog(x=flax_pops$start.long, y=flax_pops$start.lat,
                                 z=flax_pops$subpop.prevalence)
    summary(spline_prev)
    plot(spline_prev)
    # 2.5 and 97.5 percentiles represent the 95% confidence interval (very small correlation if any) 
    

# Testing for hot spots: 
    test4=lisa(x=flax_pops$start.long, y=flax_pops$start.lat, z=flax_pops$subpop.prevalence,
               neigh=5)
        # arbitrarily chose to compare 5 neighbors 
    plot(test4) 
    # would have wanted to see more filled red circles 
    
    
# Moran's test: does occurrence in one region make event in neighboring event more or less likely 
    # Adapted from https://mgimond.github.io/simple_moransI_example/
    df <- flax_pops %>% st_as_sf(coords = c("start.long", "start.lat"), dim = "XY")
    plot(df[5])
    df1 <- as(df, "Spatial")
    coo <- coordinates(df1)
    S.dist  <-  dnearneigh(coo, longlat = TRUE, 0, 1)  
        # don't actually know what constraining the boundary to be between 0 and 1 means here, 
        # just looked like it provided the most realistic neighbor connectivity when examining S.dist (~10 neighbors per pop)
    lw <- nb2listw(S.dist, style="W",zero.policy=T) 
    prev_MI  <-  moran.mc(df1$subpop.prevalence, lw, nsim=999,zero.policy=T) 
    plot(prev_MI, main="", las=1) 
    prev_MI
    # seems to be a slight relationship for prevalence,
    inc_MI  <-  moran.mc(df1$subpop.incidence, lw, nsim=999,zero.policy=T) 
    plot(inc_MI, main="", las=1) 
    inc_MI
    # No relationship for incidence 
    
    
    
    
    
# Conclusions: there might be a small relationship for prevalence but it is difficult to show statistical significance 
    # No clear relationship with incidence 
    
    


    
##  Approach 2: Using each transect as an independent survey: (currently broken, wanted to check logic with you 
    # before going further)
    
prevalence_significances <- c() 

summarized_pops<- flax_pops %>% summarize(n = length(subpop.prevalence), 
                        mu = mean(subpop.prevalence), 
                        sig = sd(subpop.prevalence) * (n - 1)/n,
                        ) 

for (transect1 in levels(flax_pops$transect)[-1]){
transect1_pop <- flax_pops %>% filter(transect == transect1)
transect1_summarized_pop <- summarized_pops %>% filter(transect == transect1)
zscale = ((transect1_pop %>% pull(subpop.prevalence)) - 
            (transect1_summarized_pop %>% pull(mu))) / 
            (transect1_summarized_pop %>% pull(sig))


# autocorrelation matrix:
rho = outer(zscale, zscale)

dst = as.matrix(dist(transect1_pop[, c("start.long", "start.lat")]))

Mantel_test <- mantel.test(M1 = rho, M2 = dst)
prevalence_significances <- c(prevalence_significances, Mantel_test$p) 

}


