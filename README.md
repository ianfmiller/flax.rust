# flax.rust
tools and analyses related to the flax-rust project

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

