# cross.scale.transmission.dynamics
This directory contains code for "Predicting the effects of climate change on the cross-scale epidemiological dynamics of a fungal plant pathogen". <br /> <br />

## data
Data files are located in <a href="https://github.com/ianfmiller/flax.rust/blob/main/data/">flax.rust/data</a>. <br> <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/pustule measurements.csv">pustule measurements.csv</a>: pustule area data <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/healthyplants.csv">healthyplants.csv</a>: uninfected focal plant data <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Withinhost.csv">Withinhost.csv</a>: infected focal plant data <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/plant locs.csv">plant locscsv</a>: plant location data <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Demography.csv">Demography.csv</a>: initial condition data for plants tracked in a parallel demographic study <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/Epidemiology.csv">Epidemiology.csv</a>: epidemiological data <br>
<a href="https://github.com/ianfmiller/flax.rust/blob/main/data/spore counts.csv">spore counts.csv</a>: spore deposition data from spore traps <br>

## analysis scripts
Analysis files are located in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics">flax.rust/cross.scale.transmission.dynamics</a>. <br><br>
The analysis pipeline begins with data cleaning. Weather data is processed in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a> , and plant location data is compiled in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant loc data set building.R">plant loc data set building.R</a>. <br /> <br />
After data cleaning, analyses begin at the within host scale. Plant growth data is compiled in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a>. We fit the GAM predicting change in plant height per day in
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a>. Using this fitted model, we generate a data set that describes the heights of all plants in each transect in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant height dataset building.R">plant height data set building.R</a> This data is used in later analyses. Next, we turn to pusutle growth. Data is prepared in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area data prep.R">pustule area data prep.R</a>. We fit the GAM predicting change in pustule area per day in
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area analysis.R">pustule area analysis.R</a>. Finally, we prepare data describing patterns of infection intensity progression in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/infection intensity data prep.R">infection intensity data prep.R</a> and fit the GAM predicting change in infection intensity per day in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/infection intensity analysis.R">infection intensity analysis.R</a>. <br /> <br />
We next analyze patterns of transmission. We fit the tilted gaussian plume model in 
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/spore deposition model fitting tilt.R">spore deposition model fitting tilt.R</a> using the functions in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/spore deposition functions.R">spore deposition functions.R</a>. We use the tilted gaussian plume model to build a data set describing predicted spore deposition, environmental conditions, and infection outcomes in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/transmission data set building.R">transmission data set building.R</a>. We fit the GAM predicting infection odds in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/transmission analysis.R">transmission analysis.R</a>.<br /> <br />
Finally, using the models fit in the above analysis scripts, we build and simulate the epidemiological model in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/simulate.epidemic.R">simulate.epidemic.R</a>

## models and summarized data

Fitted models are stored in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models">flax.rust/cross.scale.transmission.dynamics/models</a>. Processd data files are stored in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models">flax.rust/cross.scale.transmission.dynamics/summarized data</a>.
  
## figure scripts

Scripts used to generate figures are located in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/figs">flax.rust/cross.scale.transmission.dynamics/figs</a>. All simulations of trajectories from fitted GAMs occur in figure scripts.
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

