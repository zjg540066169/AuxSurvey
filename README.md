# AuxSurvey
Probability surveys often use auxiliary continuous data from administrative records, but the utility of this data is diminished when it is discretized for confidentiality. This R package provides a user-friendly function for different survey estimators with the discretized auxiliary variables.

The following survey estimators with discretized auxiliary variables are provided in the R package:
* weighted or unweighted sample mean
* weighted or unweighted raking
* weighted or unweighted poststratification
* MRP (Bayesian Multilevel regression with poststratification)
* GAMP (Bayesian Generalized additive model of response propensity)
* Bayesian linear regression
* BART (Bayesian additive regression trees)

For more details on these estimators and their applications, please consult the following paper: "Improving survey inference using administrative records without releasing individual-level continuous data".

## Installation
This package is based on `rstan` and `rstanarm`, please make sure these two packages can be installed.

Right now, the package has not been uploaded to Rcran yet, so please install this package from Github:
```
require("devtools")
install_github("https://github.com/zjg540066169/AuxSurvey")
library(AuxSurvey)
```

When you run Bayesian models, there might be the following error message. 
```
Error in initializePtr() : function 'cholmod_factor_ldetA' not provided by package 'Matrix'
```
This is because the package 'Matrix' < 1.6-2 and 'Matrix' >= 1.6-2 are binary incompatible. Run the following code to solve:
```
install.packages("lme4", type = "source")
library(AuxSurvey)
```

## Usage
There are two functions in this package. `simulate` generates datasets used in the paper. `auxsurvey` is the main function to calculate estimators. The input datasets for `auxsurvey` are data.frame or tibble. Keep all the categorical variables as factors.

### Generate simulated data
As described in paper, we generate a population dataset with 3000 samples. We then sample about 600 cases from the population dataset. These two datasets have two outcomes: continuous outcome `Y1` and binary outcome `Y2`. Covariates consist of:
* Z1 (binary): from Bernoulli(0.7).
* Z2 (binary): from Bernoulli(0.5).
* Z3 (binary): from Bernoulli(0.4).
* X (continuous): from N(0, 1).
* auX_3 (binary): discretized X with 3 categories, with equally spaced quantiles.
* auX_5 (binary): discretized X with 5 categories, with equally spaced quantiles.
* auX_10 (binary): discretized X with 10 categories, with equally spaced quantiles.
* true_pi (continuous, value in [0, 1]): the true propensity scores to generate samples.
* logit_true_pi (continuous): the logit transformation of true propensity scores.
* estimated_pi (continuous, value in [0, 1]): the estimated propensity scores. We use Bayesian additive regression trees (BART) to estimate, with 50 trees, 100 iterations for both burn-in phase and posterior sampling phase. The posterior mean is used as the estimated propensity scores.
* logit_estimated_pi (continuous): the logit transformation of estimated propensity scores.

```
library(AuxSurvey)
data = simulate(N = 3000, discretize = c(3, 5, 10), setting = 1, seed = 123) # setting parameter takes values in {1,2,3}, which corresponds the three simulation scenarios in paper.
population = data$population # get population, 3000 cases
samples = data$samples # get samples, about 600 cases
ipw = 1 / samples$true_pi # get the true inverse probability weighting in sample
est_ipw = 1 / samples$estimated_pi # get the estimated inverse probability weighting in sample
true_mean = mean(population$Y1) # true value of the estimator
```

### Estimation
After we generate datasets, we can run `auxsurvey` to get estimates. Here are the explanations of parameters:
```
auxsurvey(formula, auxiliary = NULL, samples, population = NULL, subset = NULL, family = gaussian(), method = c("sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear", "BART"), weights = NULL, levels = c(0.95, 0.8, 0.5), stan_verbose = TRUE, show_plot = TRUE, nskip = 1000, npost = 1000, nchain = 4, HPD_interval = FALSE)
```
* formula (required): A string or formula for the specified formula for the outcome model. For non-model based methods(sample mean, raking, poststratification), just include the outcome variable, such as "\~Y1" or "\~Y2". For model-based methods (MRP, GAMP, LR), additional predictors can be specified as fixed effects terms, such as "Y1\~Z1+Z2+Z3+I(Z1*Z2)". For GAMP, smooth functions can be specified, such as "Y1\~Z1+s(Z2, 10)+s(Z3, by=Z1)". Categorical variables are coded as dummy variables in model based methods.
* auxiliary (default: NULL): A string for the specified formula for the auxiliary variables. For sample mean and BART, just leave it as NULL. For raking, poststratification, GAMP, and linear regression, use string for an additive model, such as "Z1+Z2+Z3+auX_5". MRP specifies random effects for each variable (term) in this parameter, such as "Z1+Z2+Z3" or "Z1+Z2:Z3".
* samples (required): A dataframe or tibble contains all variables in 'formula' and 'auxiliary'. This dataframe is a subset of 'population'.
* population (default: NULL): A dataframe or tibble contains all variables in 'formula' and 'auxiliary'. For the sample mean estimator, it is NULL as default. If this parameter is specified for sample mean, the finite population correction will be calculate for CIs.
* subset (default: NULL): A character vector. Each element is a string representing a filtering condition to select subset of samples and population. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated. Some examples are 'c("Z1 == 1", "Z1 == 1 & Z2 == 1")'.
* family (default: `gaussian()`): The distribution family of the outcome variable. Currently, we only support `gaussian()` and `binomial()`.
* method (required): Must choose one of the following methods: "sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear", "BART".
* weights (default: NULL): Weights of each cases in `samples`.
* levels (default: c(0.95, 0.8, 0.5)): levels of CI.
* stan_verbose (default: TRUE): indicate if printing all messages when running stan models.
* show_plot	(default: TRUE): indicate if showing some diagnostic plots for stan models.
* nskip	(default: 1000): The number of burn-in iterations of each chain in MCMC for stan models.
* npost	(default: 1000): The number of posterior sampling iterations of each chain in MCMC for stan models.
* nchain (default: 4): The number of MCMC chains for stan models.
* HPD_interval (default: FALSE): indicate if the highest posterior density intervals are calculated for CI.

Here are some examples:
#### Examples: Sample mean
The parameter `population` will not be used, so just set it as NULL.
```
# Unweighted sample mean
# subset: the whole data.
sample_mean = auxsurvey("~Y1",  auxiliary = NULL, weights = NULL, samples = samples, population = NULL, subset = NULL, method = "sample_mean", levels = 0.95)

# IPW sample mean
# subset: the whole data and subset of Z1 == 1 & Z2 == 1.
# CIs are calculated with finite population correction
IPW_sample_mean = auxsurvey("~Y1",  auxiliary = NULL, weights = ipw, samples = samples, population = population, subset = c("Z1 == 1 & Z2 == 1"), method = "sample_mean", levels = 0.95)

# Estimated IPW sample mean of binary outcome
# subset: the whole data and subsets of Z1 == 1 and Z1 == 1 & Z2 == 1.
Est_IPW_sample_mean = auxsurvey("~Y2",  auxiliary = NULL, weights = est_ipw, samples = samples, population = NULL, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"), family = binomial(), method = "sample_mean", levels = 0.95)
```

#### Examples: Raking
```
# Unweighted Raking for auX_5, with interaction with Z1
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
rake_5_Z1 = auxsurvey("~Y1",  auxiliary = "Z2 + Z3 + auX_5 * Z1", weights = NULL, samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "rake", levels = 0.95)

# IPW Raking for auX_10
# subset: the whole data and subsets of Z1 == 1 and Z1 == 1 & Z2 == 1.
# CIs: 0.95, 0.8
rake_10 = auxsurvey("~Y1",  auxiliary = "Z1 + Z2 + Z3 + auX_10", weights = ipw, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"), method = "rake", levels = c(0.95, 0.8))

# Estimated IPW Raking for auX_3, binary outcome
# subset: the whole data.
rake_3 = auxsurvey("~Y2",  auxiliary = "Z1 + Z2 + Z3 + auX_3", weights = est_ipw, samples = samples, population = population, subset = NULL, family = binomial(), method = "rake", levels = 0.95)
```

#### Examples: Poststratification
Since auX_5 and auX_10 will lead to 40 and 80 stratas, which are too many. We only applied to auX_3.
```
# Unweighted Poststratification for auX_3
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
ps_3 = auxsurvey("~Y1",  auxiliary = "Z1 + Z2 + Z3 + auX_3", weights = NULL, samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "postStratify", levels = 0.95)

# IPW Poststratification for auX_3, binary outcome
# subset: the whole data and subsets of Z1 == 1 and Z1 == 1 & Z2 == 1.
# CIs: 0.95, 0.8
ps_3_binary = auxsurvey("~Y2",  auxiliary = "Z1 + Z2 + Z3 + auX_3", weights = ipw, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"), family = binomial(), method = "postStratify", levels = c(0.95, 0.8))
```

#### Examples: MRP
MRP treats `formula` as the fixed effects part, and treats `auxiliary` as the random effects part. `weights` is not used for MRP. MRP allows the nested models.
```
# MRP with auX_3
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
# The model is Y1 ~ 1 + Z1 + Z2 + (1|Z3) + (1|auX_3)
MRP_1 = auxsurvey("Y1~1 + Z1 + Z2",  auxiliary = "Z3 + auX_3", samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "MRP", levels = 0.95)

# MRP with auX_10, nested within Z1
# subset: the whole data.
# The model is Y1 ~ 1 + Z1 + Z2 + Z3 + (1|Z1:auX_3)
MRP_2 = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "Z1:auX_3", samples = samples, population = population, subset = NULL, method = "MRP", levels = 0.95)

# MRP with auX_10, nested within Z1. Z3 will be used in both fixed and random parts. Outcome is binary.
# subset: the whole data.
# CI are HPD_interval
# The model is Y2 ~ 1 + Z1 + Z2 + Z3 + (1|auX_3) + (1|Z3).
MRP_3 = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxiliary = "Z3 + auX_3", samples = samples, population = population, subset = NULL, family = binomial(), method = "MRP", levels = 0.95, HPD_interval = T)
```
#### Examples: GAMP
GAMP specifies smooth functions for discretized variables and other variables. Smooth functions can be used in both `formula` and `auxiliary` parameters. If variables are specified in `auxiliary` without smooth function, random effects will be used.
```
# GAMP with smooth functions on auX_3 and logit_estimated_pi
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
# The model is Y1 ~ 1 + Z1 + Z2 + Z3 + s(logit_estimated_pi) + s(auX_3)
GAMP_1 = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "s(logit_estimated_pi) + s(auX_3)", samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "GAMP", levels = 0.95)
# This is equivalent to
# GAMP_1 = auxsurvey("Y1~1 + Z1 + Z2 + Z3 + s(logit_estimated_pi)",  auxiliary = "s(auX_3)", samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "GAMP", levels = 0.95)


# GAMP with smooth functions on logit_estimated_pi and auX_10 with interaction on Z1
# subset: the whole data.
# The model is Y1 ~ 1 + Z1 + Z2 + Z3 + s(logit_estimated_pi, by = Z1) + s(auX_3, by = Z1)
GAMP_2 = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "s(logit_estimated_pi, by = Z1) + s(auX_3, by = Z1)", samples = samples, population = population, subset = NULL, method = "GAMP", levels = 0.95)


# GAMP with smooth functions on logit_estimated_pi and auX_10 with interaction on Z1, and Z3 as random effects
# subset: the whole data.
# CI are HPD_interval
# The model is Y1 ~ 1 + Z1 + Z2 + s(logit_estimated_pi, by = Z1) + s(auX_3, by = Z1) + (1|Z3)
GAMP_3 = auxsurvey("Y1~1 + Z1 + Z2",  auxiliary = "s(logit_estimated_pi, by = Z1) + s(auX_3, by = Z1) + Z3", samples = samples, population = population, subset = NULL, method = "GAMP", levels = 0.95, HPD_interval = T)

# GAMP with smooth functions on logit_estimated_pi and auX_10 with interaction on Z1, and Z1:Z3 as random effects. We can specify the number of degrees of smooth function. Outcome is binary.
# subset: the whole data.
# The model is Y2 ~ 1 + Z1 + Z2 + s(logit_estimated_pi, by = Z1, k = 3) + s(auX_3, by = Z1, k = 5) + (1|Z1:Z3)
GAMP_4 = auxsurvey("Y2~1 + Z1 + Z2",  auxiliary = "s(logit_estimated_pi, by = Z1, k = 3) + s(auX_3, by = Z1, k = 2) + Z1:Z3", samples = samples, population = population, subset = NULL, family = binomial(), method = "GAMP", levels = 0.95)
```

#### Examples: Bayesian Linear Regression
Bayesian linear regression treats all categorical variables as dummy variables. `auxiliary` will be directly added to `formula`.
```
# Linear regression
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
# CI are HPD_interval
# The model is Y1 ~ 1 + Z1 + Z2 + Z3 + auX_3.
LR_1 = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "auX_3", samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "linear", levels = 0.95,  HPD_interval = T)

# Linear regression with quadratic terms of auX_5
# subset: the whole data.
# The model is Y1 ~ 1 + Z1 + Z2 + Z3 + auX_5 + auX_5^2
LR_2 = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "auX_5 + I(auX_5^2)", samples = samples, population = population, subset = NULL, method = "linear", levels = 0.95)
```

#### Examples: Bayesian additive regression trees
For BART, specify the predictors in `formula` and leave `auxiliary` as NULL.
```
# BART
# subset: the whole data and subset of Z1 == 1 & Z2 == 1 & Z3 == 0.
# CI are HPD_interval
BART_1 = auxsurvey("Y1~ Z1 + Z2 + Z3 + auX_3",  auxiliary = NULL, samples = samples, population = population, subset = "Z1 == 1 & Z2 == 1 & Z3 == 0", method = "BART", levels = 0.95,  HPD_interval = T)
```



## Disclaimer

If you find there is any bug, please contact me: jungang.zou@gmail.com.
