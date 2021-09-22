# cross.scale.transmission.dynamics
This file contains analyses for a paper (in preparation) about the effects of climate change on the cross-scale transmission dynamics of flax rust disease. <br /> <br />
The below text provides a summary of the analysis pipeline. To recreate results source the following files in order: <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a> <br />
<a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a> <br />
## Data collection
## Statistical methods
Unless otherwise noted, we use the following modeling approach to investigate the relationship between a response variable (e.g. observed plant height at time t+1) and predictor variables. These predictor variables include the previous observed state (e.g. plant height at time t) and weather metrics (e.g. mean temperature).
<br />
<br />
We use generalized additive models (implemented via the mgcv package) to model the response variable as the sum of smooth functions of predictors: y ~ s(x<sub>1</sub>) + s(x<sub>2</sub>) + ... + s(x<sub>n</sub>). These smooths are capable of capturing non-linear relationships between the predictors and the response variable. We limit the number of knots in each smooth to four in order to prevent overfitting. We include interactions with elapsed time for all predictor variables (except total rainfall) to account for the variation in the elapsed time between observations. To select the most parsimonious model we begin by fitting a model with all possible predictors using penalized cubic regression splines with an additional shrinkage penalty. This approach functionally removes insignificant predictors from the model, facilitating predictor selection (see <a href="https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html">mgcv documentation</a>). Next, we iteratively remove insignificant terms from the model until none remained.
## Weather data prep
Weather data collected from sensors at each of four sites is compiled and cleaned in <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/prep.enviro.data.R">prep enviro data.R</a>. 
## Within host scale
The analysis begins at the within host scale. 
### Plant growth
First, we investigate the relationship between plant growth and weather conditions.
<br />
<br />
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth data prep.R">plant growth data prep.R</a> longitudinal height data of healthy and diseased focal plants is joined. This raw data is stored as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS">plant.heights.RDS</a>. To make this data usable for analyses, we join data on change in plant height with mean weather metrics in a new data object <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/summarized data/delta.heights.RDS">delta.heights.RDS</a>. 
<br />
<br />
In <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/plant growth analysis.R">plant growth analysis.R</a> we fit and visualize a model of plant height<sub>t+1</sub> as a sum of smoothed functions of plant height<sub>t</sub> and weather metrics (interactions with time are included for all predictors). This model is saved as <a href="https://github.com/ianfmiller/flax.rust/blob/main/cross.scale.transmission.dynamics/models/plant.growth.model.RDS">plant.growth.model.RDS</a>.
### Pustlule growth
-pustule area data prep.R <br />
-pustule area analysis.R <br />
-within host climate prediction functions.R (vis)
### Pustule establishment
-n pustules data prep.R <br />
-n pustules analysis.R <br />
-within host climate prediction functions.R (vis)
### Plant infection intensity change
-plant data prep.R <br />
-plant analysis.R <br />
-predict plant inf intens change funcs.R <br />
-within host climate prediction functions.R (vis)
## Between host scale
### Data prep
#### Plant heights
-plant height dataset building.R
#### Plant locations and incidence data
-plant loc dataset building.R
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

