# AuxSurvey
Probability surveys often use auxiliary continuous data from administrative records, but the utility of this data is diminished when it is discretized for confidentiality. This R package provides a user-friendly function for different survey estimators with the discretized auxiliary variables.

The following survey estimators with discretized auxiliary variables are provided in the R package:
* weighted or unweighted sample mean
* weighted or unweighted raking
* weighted or unweighted poststratification
* MRP (Bayesian Multilevel regression with poststratification)
* GAMP (Bayesian Generalized additive model of response propensity)
* Bayesian linear regression

For more details on these estimators and their applications, please consult the following paper: "Improving survey inference using administrative records without releasing individual-level continuous data".

## Installation
This package is based on `rstan` and `rstanarm`, please make sure these two packages can be installed.

Right now, the package has not been uploaded to Rcran yet, so please install this package from Github:
```
require("devtools")
install_github("https://github.com/zjg540066169/AuxSurvey")
```

## Usage
There are two functions in this package. `simulate` generates datasets used in the paper. `auxsurvey` is the main function to calculate estimators.

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
```

### Estimation
After we generate datasets, we can run `auxsurvey` to get estimates. Here are the explanations of parameters:
```
auxsurvey(formula, auxiliary = NULL, samples, population = NULL, subset = NULL, family = gaussian(), method = c("sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear"), weights = NULL, levels = c(0.95, 0.8, 0.5), stan_verbose = TRUE, show_plot = TRUE, nskip = 1000, npost = 1000, nchain = 4, HPD_interval = FALSE)
```
* formula (required): A string or formula for the specified formula for the outcome model. For non-model based methods(sample mean, raking, poststratification), just include the outcome variable, such as "\~Y1" or "\~Y2". For model-based methods (MRP, GAMP, LR), additional predictors can be specified as fixed effects terms, such as "Y1\~Z1+Z2+Z3+I(Z1*Z2)". For GAMP, smooth functions can be specified, such as "Y1\~Z1+s(Z2, 10)+s(Z3, by=Z1)". Categorical variables are coded as dummy variables in model based methods.
* auxiliary (default: NULL): A string for the specified formula for the auxiliary variables. For sample mean, just leave it as NULL. For raking, poststratification, GAMP, and linear regression, use string for an additive model, such as "Z1+Z2+Z3+auX_5". MRP specifies random effects for each variable (term) in this parameter, such as "Z1+Z2+Z3" or "Z1+Z2:Z3".
* samples (required): A dataframe or tibble contains all variables in 'formula' and 'auxiliary'. This dataframe is a subset of 'population'.
* population (default: NULL): A dataframe or tibble contains all variables in 'formula' and 'auxiliary'. For the sample mean estimator, it doesn`t need information from population, so this parameter is NULL in this case.
* subset (default: NULL): A character vector. Each element is a string representing a filtering condition to select subset of samples and population. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated. Some examples are 'c("Z1 == 1", "Z1 == 1 & Z2 == 1")'.
* family (default: `gaussian()`): The distribution family of the outcome variable. Currently, we only support `gaussian()` and `binomial()`.
* method (required): Must choose one of the following methods: "sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear".
* weights (default: NULL): Weights of each cases in `samples`.
* levels (default: c(0.95, 0.8, 0.5)): levels of CI.
* stan_verbose (default: TRUE): indicate if printing all messages when running stan models.
* show_plot	(default: TRUE): indicate if showing some diagnostic plots for stan models.
* nskip	(default: 1000): The number of burn-in iterations of each chain in MCMC for stan models.
* npost	(default: 1000): The number of posterior sampling iterations of each chain in MCMC for stan models.
* nchain (default: 4): The number of MCMC chains for stan models.
* HPD_interval (default: FALSE): indicate if the highest posterior density intervals are calculated for CI.

```
# shrinkage models
model1 = Laplace(Y_array, X_array, standardize = True, r = 2, s = 15)
# model1 = Horseshoe(Y_array, X_array, standardize = True)
# model1 = Ridge(Y_array, X_array, standardize = True)

# discrete mixture models
model2 = Spike_laplace(Y_array, X_array, standardize = True, lambda_ = 6/11, a = 1, b = 1)
# model2 = Spike_ridge(Y_array, X_array, standardize = True, p0 = 0.5, v0 = 4)
```
Here `Y_array` is a 2-d data array for response variable, its dimension is `(n_imputations, n_samples)`. `X_array` is a 3-d data array for explanatory variables, its dimension is `(n_imputations, n_samples, n_features)`. If the parameter `standardize` is True, X_array will be standardized and then used to run MCMC chains. If it is False, the original X_array is used to calculate MCMC chains. Other parameters are hyper-parameters for each model.


### Posterior Sampling
After initialization, we can use `sample` function to run MCMC chains and get posterior samples:
```
model1.sample(n_post = 1000, n_burn = 500, target_accept = 0.9, n_chain = 2, n_thread = 4, max_treedepth = 10, seed = 123)
# model2.sample(n_post = 1000, n_burn = 500, target_accept = 0.9, n_chain = 2, n_thread = 4, max_treedepth = 10, seed = 123)
```
The parameters for `sample` function are as follows:
* n_post(required): number of posterior samples for each chain.
* n_burn(required): number of burn-in samples for each chain.
* target_accept(default 0.9): target acceptance probability for NUTS.
* max_treedepth(default 10): maximum tree depth for NUTS.
* n_chain(default 1): number of parallel chains to run.
* n_thread(default 4): number of threads to run parallel computing.
* seed(default None): random seed. If seed is None, seed is equals to the current time in seconds since the Epoch.

We can use `get_posterior_samples` function to get posterior samples:
```
model1.get_posterior_samples(var_name = "beta", rescale = True)   # get posterior samples for the coefficients vector
model2.get_posterior_samples(var_name = "alpha", rescale = True)  # get posterior samples for the intercept
model2.get_posterior_samples(var_name = "g", rescale = True)      # get posterior samples for the hidden binary variables in discrete mixture models
```
Here `var_name` is the variable we want to sample for. `rescale` specifies whether to return coefficients in the original scale; if it is False, then coefficients corresponding to standardized covariates are return; if it is True, all the coefficients are rescaled to the original scale. If `standardize = False` in initialization stage, `rescale` has no effect. For MI data, we simply mixed up the posterior samples for each grouped coefficient among all MI sets. So the dimension of posterior samples for coefficients vector `beta` is `(n_chains, n_imputations * n_samples, n_features)`. And the dimension of intercept `alpha` is `(n_chains, n_imputations * n_samples)`.

### Summary Statistics
Our library provides a `summary` function to generate summary statistics for all the variables in the hierachical model:
```
summary_stats1 = model1.summary(rescale = True)
print(summary_stats1)
```
Here `rescale` is the same as it in function `get_posterior_samples`.


### Variable Selection
Users can use `select` function to select important variables:
```
select1 = model1.select(value = 0.95, rescale = True) # Credible Interval Criterion for Shrinkage Models
select2 = model2.select(value = 0.5,  rescale = True) # Cutting-off point for Discrete Mixture Models
```
The meaning of `value` depends on the type of models. For shrinkage models, `value` is the credible interval criterion for selection. For discrete mixture models, `value` stands for the cutting-off point for selection. For more details, please consult Chapter 3.2 in the paper: "Variable Selection for Multiply-imputed Data: A Bayesian Framework" (Arxiv: https://arxiv.org/abs/2211.00114).


### Evaluation
There are some evaluation functions in this library:
```
from bmiselect.utils.evaluation import sensitivity, specificity, f1_score, mse

sensitivity(select = select1, truth = truth)                                           # sensitivity
specificity(select = select2, truth = truth)                                           # specificity
f1_score(select = select2, truth = truth)                                              # f1 score
mse(beta, covariance, select = select1, X = X_array, Y = Y_array, intercept = True)    # mse, given coefficients and covariance matrix of ground truth
```
Here `select` and `truth` are binary vectors with length `(n_features)`. `select[i] = True` means i-th variable is selected.


### Refitting Linear Regression
After we complete the variable selection by Bayesian MI-LASSO, users can apply `fit_lr` to fit ordinary linear regression separately on each imputed dataset. Alternatively, users can utilize `pooled_coefficients` and `pooled_covariance` to directly get pooled coefficients and covariance matrix with Rubin`s Rule.
```
from bmiselect.utils.evaluation import fit_lr, pooled_coefficients, pooled_covariance

# refit linear regression by using selected variabels
lr_models = fit_lr(select1, X_array, Y_array, intercept = True))
for lr in lr_models:
    print(lr.summary())

# get pooled coefficients estimates by using Rubin`s rule
lr_coef = pooled_coefficients(select2, X_array, Y_array, intercept = True))
print(lr_coef)

# get pooled covariance estimates by using Rubin`s rule
lr_covariance = pooled_covariance(select2, X_array, Y_array, intercept = True))
print(lr_covariance)
```
If `intercept = True`, then an intercept is added to each ordinary linear regression respecitively. If `intercept = False`, then no intercept is used in linear regression.


## Disclaimer

If you find there is any bug, please contact me: jungang.zou@gmail.com.
