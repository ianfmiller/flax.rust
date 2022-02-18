# cross.scale.transmission.dynamics
This file contains code for "Predicting the effects of climate change on the cross-scale epidemiological dynamics of a fungal plant pathogen". <br /> <br />

## data
Data files are located in the <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/">data directory.R</a>. Pustule area data is given in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/pustule measurements.csv">pustule measurements.csv</a>. Data from uninfected focal plants is contained in in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/healthyplants.csv">healthyplants.csv</a>, and data from infectd focal plants is contained in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Withinhost.csv">Withinhost.csv</a>. Plant location data is given in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/plant locs.csv">plant locscsv</a>. <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Demography.csv">Demography.csv</a> describes the first observed conditions of several plants tracked as part of a parallel demographic study. <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Epidemiology.csv">Epidemiology.csv</a> conntains data on the spatiotemporal spread of infection. Spore deposition data from spore traps is given in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/spore counts.csv">spore counts.csv</a>.

## analysis scripts
The analysis pipeline begins with data cleaning. Weather data is processed in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a> , and plant location data is compiled in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant loc data set building.R">plant loc data set building.R</a>. <br /> <br />
After data cleaning, analyses begin at the within host scale. Plant growth data is compild in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a> and analyzed in
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a>. Using the fitted GAM describing plant growth, we generate a data set in that describes the heights of all plants in each transect. <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant height dataset building.R">plant height data set building.R</a>  This data is used in later analyses. Next, we turn to pusutle growth. Data is prepared in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area data prep.R">pustule area data prep.R</a> and analyzed in
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area analysis.R">pustule area analysis.R</a>. Finally, we prepare data describing patterns of infection intensity progression in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/infection intensity data prep.R">infection intensity data prep.R</a> and analyze this data in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/infection intensity analysis.R">infection intensity analysis.R</a>. <br /> <br />
We next analyze patterns of transmission. We fit the tilted gaussian plume model in 
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/spore deposition model fitting tilt.R">spore deposition model fitting tilt.R</a> using the functions in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/spore deposition functions.R">spore deposition functions.R</a>. We use the tilted gaussian plume model to build a data set describing predicted spore deposition, environmental conditoins, and infection outcomes in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/transmission data set building.R">transmission data set building.R</a>. We analyze this data in  <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/transmission analysis.R">transmission analysis.R</a>.<br /> <br />
Finally, using the models fit in the above analysis scripts, we build and simulate the epidemiological model in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/simulate.epidemic.R">simulate.epidemic.R</a>

## models and summarized data

## figure scripts

# quant.vir.evolution

## infeciton intensity survival analysis--demog data only.R
### Uses only 1st observations from within a season as predictors
### predicting probability of death from disease metrics and size metrics
### conclusions:
% stems infected increases odds of death; increased height decreases odds of death. <br />
N stems infected increases odds of death; total N stems decreases odds of death.

## infeciton intensity survival analysis
### predict death using only 1st observations of plants within a season
part 1: disease metrics as predictors, including healhty plant data <br />
part 2: disease metrics and size metrics as predictors, including H+D plant data <br />
part 3: disease metrics and size metrics at 1st observation as predictors, incluuding H+D plant data
#### conclusions:
Infection slightly increases odds of death for larger plants

## survival analysis
### Uses only 1st observations from within a season as predictors
### predicting probability of death from disease status and size metrics
### conclusions:
N stems decreases odds of death; no significant effect of disease status

## categorical outcome analysis
### ordinal logistic regression to predict status (H,D,X) from disease/size metrics
### conclusions: more diseased, small size correlated with D, then X, but effects are just barely statistically insignificant

