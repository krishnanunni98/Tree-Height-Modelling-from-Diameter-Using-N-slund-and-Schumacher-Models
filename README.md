# Tree Height Modelling from Diameter Using Näslund and Schumacher Models

## Project Overview

This exercise models tree height as a function of diameter at breast height using classical Näslund height–diameter equations and a Schumacher-type mixed-effects model. The workflow includes data preparation, plotwise OLS and nonlinear estimation, bias-corrected prediction, and graphical comparison of fitted curves. 

## Objectives

- Model the tree height–diameter relationship.
- Compare linearised and nonlinear Näslund models.
- Fit a mixed-effects Schumacher model with plot-level random effects.
- Apply bias correction and back-transformed prediction.
- Visualize fitted curves against observed tree data. 

## Data

- `tree_list.dat`
- Variables: `plot_id`, `tree_id`, `plot.area.m2`, `T.yrs`, `d13.cm`, `h.m` :contentReference[oaicite:2]{index=2}

## Methodology

### Data Preparation
The tree list was ordered by plot and diameter, sample trees were separated, and transformed predictors were created for height modelling. :contentReference[oaicite:3]{index=3}

### Näslund Modelling
A linearised Näslund model was fitted with OLS for a selected plot, and then a nonlinear Näslund model was fitted with `nls()`. The exercise diary shows both the linear and nonlinear fitted curves together with the observed height–diameter data. 

### Mixed-Effects Modelling
A Schumacher-style model was fitted with `nlme::lme()`, using `plot_id` as a random effect to account for between-plot variation. The model was then used to produce localized, bias-corrected height predictions. :contentReference[oaicite:5]{index=5}

## Main Outputs

- Height–diameter scatter plots
- Näslund linear and nonlinear fitted curves
- Schumacher mixed-effects height curve
- Bias-corrected localised predictions
- Parameter estimates and variance components 

## Skills Demonstrated

- Tree height modelling
- Height–diameter curve fitting
- Linear regression and nonlinear regression
- Mixed-effects modelling
- Bias correction and back-transformation logic
- Plot-level random effects
- Interpretation of inventory model outputs 

