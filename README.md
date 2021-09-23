# cross.scale.transmission.dynamics
This file contains analyses for a paper (in preparation) about the effects of climate change on the cross-scale transmission dynamics of flax rust disease. <br /> <br />
The below text provides a summary of the analysis pipeline. To recreate results, source the following files in order: <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant loc data set building.R">plant loc data set building.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant height dataset building.R">plant height data set building.R</a>  <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area data prep.R">pustule area data prep.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area analysis.R">pustule area analysis.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/n pustules data prep.R">n pusutles data prep.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/n pustules analysis.R">n pustules analysis.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant data prep">plant data prep.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant analysis.R">plant analysis.R</a> <br />
## Data collection
## Statistical methods
Unless otherwise noted, we use the following modeling approach to investigate the relationship between a response variable (e.g. observed plant height at time t+1) and predictor variables. These predictor variables include the previous observed state (e.g. plant height at time t) and weather metrics (e.g. mean temperature).
<br />
<br />
We use generalized additive models (implemented via the mgcv package) to model the response variable as the sum of smooth functions of predictors: y ~ s(x<sub>1</sub>) + s(x<sub>2</sub>) + ... + s(x<sub>n</sub>). These smooths are capable of capturing non-linear relationships between the predictors and the response variable. We limit the number of knots in each smooth to four in order to prevent overfitting. We include interactions with elapsed time for all predictor variables (except total rainfall) to account for the variation in the elapsed time between observations. To select the most parsimonious model, we begin by fitting a model with all possible predictors using penalized cubic regression splines with an additional shrinkage penalty. This approach functionally removes insignificant predictors from the model, facilitating predictor selection (see <a href="https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html">mgcv documentation</a>). Next, we iteratively remove insignificant terms from the model until none remained.
## General data prep
### Weather data
Weather data collected from sensors at each of four sites is compiled and cleaned in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a>. 
### Plant locations and epidemiological data
We recorded the locations of each plant in each transect at the begininning of the field season. When we observed a newly diseased plant, we re-recorded its location. In order to understand spatial epidemiological patterns, we need to match newly diseased plants to the healthy plants surveyed at the beginning of the season. This is accomplished in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant loc data set building.R">plant loc data set building
  
  
  
  
  
  
  
  
  
  
  
  
  
  .R</a>. This script generates two data sets. The first, <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/corrected.locs.RDS">corrected.locs.RDS</a>, gives the identity location of all (healthy and diseased) plants in each transect. The second, <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/corrected.locs.RDS">corrected.locs.RDS</a>, gives the identity, location, and date first observed diseased for all infected plants. When the newly diseased plant was previously tagged, we updated its coordinates in the epidemiological data set to match its coordinates in the plant location data set. When the newly diseased plant was not previously tagged (and therefore able to be definitively matched to a record in the initial transect survey) we matched its identity and changed its coordinates to that of the closest untagged and unmatched healthy plant in the plant location daataset. Note that this procedure was also performed for plants that were tagged after being observed diseased. When no untagged and unmatched plant was withing 0.5m of the newly diseased plant we assumed that it was either unemerged or missed during the initial survey, and added a new record into the plant location dataset to reflect its presence.
-plant loc dataset building.R
## Within host scale
The analysis begins at the with in host scale. 
### Plant growth
First, we investigate the relationship between plant growth and weather conditions.
#### Data preperation
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a> longitudinal height data of healthy and diseased focal plants is joined. This raw data is stored as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS">plant.heights.RDS</a>. To make this data usable for analyses, we join data on change in plant height with mean weather metrics in a new data object <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/delta.heights.RDS">delta.heights.RDS</a>. As a part of this data preparation, raw data on plant infection intensity is also compiled as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/plants.RDS">plants.RDS</a>
#### Model fitting
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a> we fit and visualize a model of plant height<sub>t+1</sub> as a sum of smoothed functions of plant height<sub>t</sub> and weather metrics (interactions with time are included for all predictors), and proceeded with our standard model selection approach. This model is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/plant.growth.model.RDS">plant.growth.model.RDS</a>.
#### Plant height dataset building
In a subsequent analysis of the effects of spore deposition on the odds of infection, we wish to use plant height as a predictor variable. However, for the majority of plants within the transect, heights were not recorded after the initial population survey. We use the <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/plant.growth.model.RDS">plant.growth.model.RDS</a> to simulate missing height data from weather conditions and the closest recorded height measurment for each plant in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant height dataset building.R">plant height data set building.R</a>. This script uses convenience functions in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/predict plant height change funcs.R">predict plant height change funcs.R</a>. The resutling dataset of observed and fore/hindcasted plant heights is stored as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/corrected.plant.heights.RDS">corrected.plant.heights.RDS</a>.
### Pustlule growth
Next, we begin our analysis of the effects of climate change on the within host spread of disease beginning that the scale of a single pustule. 
#### Data preperation
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area data prep.R">pustule area data prep.R</a> we calculate pustule areas from the minimum and maximum diameters measured from pictures taken in the field. This data, along with other metadata including the date and time at which the pictures were taken are saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/pustules.RDS">pustules.RDS</a>. To facilitate analyses, we reformat this data, joining data on change in pustule diameter with mean weather metrics. This new data object is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS">delta.pustules.RDS</a>. 
#### Model fitting
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/pustule area analysis.R">pustule area analysis.R</a> we fit a model of pustule area<sub>t+1</sub> as a sum of smoothed functions of pustule area<sub>t</sub> and weather metrics (interactions with time are included for all predictors), and proceeded with our standard model selection approach. This model is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/pustule.model.RDS">pustules.model.RDS</a>. 
#### Model interpretation and visualization
In order to translate the fitted model into predictions about climate effects on pustule growth, we need to be able to make predictions from the model. This requires constructing a dummy data set for the starting pustule area and weather conditions that we wish to investigate. These dummy data sets are constructed using convenience functions in  <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/within host climate prediction functions">within host climate prediction functions</a>. Using these dummy datasets, we make predictions via bootstrap simulations of the fitted model. 
<br />
<br />
We take two approaches to visualize the effect of climate change on the growth of pustules of various size. In the first, we simulate pustule growth using the fitted model and real weather data from days corresponding to the 50th, 75th, and 90th quantiles of hottest days observed during the data collection period. In the second, we simulate pustule growth using the fitted model and weather data from a day corresponding to 50th quantile of hottest days observed during the data collection period with either 0, 1.8, or 3.7 degrees added to individual temperature observations. In each approach we simulate area change across a one day window. 
<br />
<br />
To predict how climate change will affect the growth trajectory of pustules, we use the fitted model and longitudinal weather data to simulate the growth of a pustule over time. We added either 0, 1.8 or 3.7 degrees to observed temperature readings to investigate the effects of climate change. We assumed a starting pustule area of .01cm<sup>2</sup>, and simulated trajectories for 100 pustules at each of the temperature scenarios. We replicated this procedure for forward simulation windows of one, two, and seven days. 
### Pustule establishment
We next consider the process of pustule establishment. The main metric we use to investigate this process is the change in the number of pustules present on a leaf.
#### Data preperation
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/n pustules data prep.R">n pusutles data prep.R</a> we clean raw pustule count data, and save it as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS">n.pustules.RDS</a>. To facilitate analyses, we reformat this data, joining data on change in pustule number with mean weather metrics and a predictor variable capturing how the observed weather conditions are expected to affect pustule growth. To generate this additional variable, we predicted pustule area<sub>t+1</sub> using the <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/pustule.model.RDS">pustules.model.RDS</a> model and weather metrics spanning t to t+1, while assuming pustule area<sub>t</sub>=0.01 cm<sup>2</sup>.  The new data object is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS">delta.n.pustules.RDS</a>. 
#### Model fitting
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/n pustule analysis.R">n pustules analysis.R</a> we fit a model of N pustules<sub>t+1</sub> as a sum of smoothed functions of N pustules<sub>t</sub>, predicted pustule growth, and weather metrics (interactions with time are included for all predictors), and proceeded with our standard model selection approach. This model is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/pustule.model.RDS">pustules.model.RDS</a>. 
#### Model interpretation and visualization
We visualize the effects of climate change on pustule establishment by simulating <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/pustule.model.RDS">pustules.model.RDS</a>, while again making use of convenience functions in   <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/within host climate prediction functions">within host climate prediction functions</a>. First, we simulated how climate change would affect pustule establishment for leaves with various numbers of pustules already present. We considered the same two sets of climate scenarios descrbied in the pustule area section. Next, we use the fitted model and longitudinal weather data to simulate the change in the number of pustules on one leaf over time using the same procedure outlined in the pustule area section We assumed that leaves started with one pustule. 
### Plant infection intensity change
The final within host scale we consider is the intensity of infection within an entire plant. We construct a metric to measure this intensity by multiplying the number of diseased stems by the average length of diseased tissue on a stem by the average number of pustules found on a leaf in the middle of the region of infected tissue. Note that this is not an estimate of the total number of pustules on a plant as it ignores pustules present on stems.
#### Data preperation
Raw data on plant infection intensity is summarized into <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/plants.RDS">plants.RDS</a> in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a>. This data is reformated in  <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant data prep">plant data prep.R</a> so that infection intensities at time t+1 are paired with infection intensities at time t and summary weather metrics. This new data object is stored as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS">delta.plants.RDS</a>.
#### Model fitting
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant analysis.R">plant analysis.R</a> we fit a model of infection intensity<sub>t+1</sub> as a sum of smoothed functions of infection intensity<sub>t</sub> and weather metrics (interactions with time are included for all predictors). We included 10 knots for the smoothed function of infection intensity<sub>t</sub> in order to account for the signifincant nonlineaar relationship between this term and infection intensity<sub>t+1</sub>. Our model selection approach followed the standard procedure outlined above. This model is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/plants.model.RDS">plants.model.RDS</a>. 
#### Model interpretation and visualization
We visualize the effects of climate change on within plant infection intensity by simulating <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/pustule.model.RDS">plants.model.RDS</a>, while again making use of convenience functions in   <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/within host climate prediction functions">within host climate prediction functions</a>. First, we simulated how climate change would affect the change in infection intnsity for plants with infections of varying intensity. We considered the same two sets of climate scenarios descrbied in the pustule area section. Next, we use the fitted model and longitudinal weather data to simulate the change in plant infection intensity over time using the same procedure outlined in the pustule area section. We assumed plants to begin with an infection intensity score of 10.
## Between host scale
### Spore deposition
-spore deposition functions tilt.R <br />
-spore deposition model fitting tilt.R <br />
-visualize.spore.kernel.R (vis)
### Spore deposition and infection
-building foi dataset.R <br />
-foi analysis.R
### Epidemic simulation
-epi data vis.R <br />
-starting plant inf intens model.R <br />
-simulate epidemic.R

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

