# flax.rust
tools and analyses related to the flax-rust project

## cross.scale.transmission.dynamics

### Weather data prep
-prep.enviro.data.R
### Within host scale
#### Plant growth
-plant growth data prep.R\n
-plant.growth.analysis.R
#### Pustlule growth
-pustule area data prep.R\n
-pustule area analysis.R\n
-within host climate prediction functions.R (vis)
#### Pustule establishment
-n pustules data prep.R
-n pustules analysis.R
-within host climate prediction functions.R (vis)
#### Plant infection intensity change
-plant data prep.R
-plant analysis.R
-predict plant inf intens change funcs.R
-within host climate prediction functions.R (vis)
### Between host scale
#### Data prep
##### Plant heights
-plant height dataset building.R
##### Plant locations and incidence data
-plant loc dataset building.R
#### Spore deposition
-spore deposition functions tilt.R
-spore deposition model fitting tilt.R
-visualize.spore.kernel.R (vis)
#### Spore deposition and infection
-building foi dataset.R
-foi analysis.R
#### Epidemic simulation
-epi data vis.R
-starting plant inf intens model.R
-simulate epidemic.R

## quant.vir.evolution

### infeciton intensity survival analysis--demog data only.R
#### Uses only 1st observations from within a season as predictors
#### predicting probability of death from disease metrics and size metrics
#### conclusions:
% stems infected increases odds of death; increased height decreases odds of death.
N stems infected increases odds of death; total N stems decreases odds of death.

### infeciton intensity survival analysis
#### predict death using only 1st observations of plants within a season
part 1: disease metrics as predictors, including healhty plant data
part 2: disease metrics and size metrics as predictors, including H+D plant data
part 3: disease metrics and size metrics at 1st observation as predictors, incluuding H+D plant data
#### conclusions:
Infection slightly increases odds of death for larger plants

### survival analysis
#### Uses only 1st observations from within a season as predictors
#### predicting probability of death from disease status and size metrics
#### conclusions:
N stems decreases odds of death; no significant effect of disease status

### categorical outcome analysis
#### ordinal logistic regression to predict status (H,D,X) from disease/size metrics
#### conclusions: more diseased, small size correlated with D, then X, but effects are just barely statistically insignificant

