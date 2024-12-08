library(mgcv)
library(rstanarm)

#' Simulate Survey Data with Discretized Auxiliary Variables
#'
#' This function simulates survey data with discretized auxiliary variables. It generates a population
#' dataset with continuous and binary outcomes, and includes auxiliary variables that are discretized
#' into multiple categories. The function also generates a subset of the population as a sample, based on
#' the propensity scores.
#'
#' The function supports multiple simulation settings, where each setting modifies the relationships
#' between the outcome variables and the covariates.
#'
#' @param N Number of population units to simulate. Default is 3000.
#' @param discretize A scale specifying the number of categories for discretizing continuous variables.
#'                   The function discretizes both `X` and `W` into the specified categories. Default is a number among `(3, 5, 10)`.
#' @param setting A numeric value to specify the simulation setting. The settings define different relationships
#'                between the outcome variables and the covariates. Possible values are 1, 2, 3, and 4.
#'                Default is a number among `c(1, 2, 3)`.
#' @param seed An optional random seed for reproducibility. Default is `NULL`.
#'
#' @return A list containing two elements:
#' - `population`: A tibble with the simulated population data, including both continuous and binary outcomes,
#'   as well as auxiliary variables (both raw and discretized).
#' - `samples`: A tibble with the simulated sample data, where individuals are included based on their estimated
#'   propensity scores.
#'
#' @examples
#' # Simulate survey data with setting 1 and discretizing variables 3 categories
#' data = simulate(N = 3000, discretize = 3, setting = 1, seed = 123)
#'
#' # Extract population and sample datasets
#' population = data$population
#' samples = data$samples
#'
#' # Examine the simulated population data
#' head(population)
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom gtools quantcut
#' @import mgcv
#' @importFrom dplyr filter %>% as_tibble
#' @import stringr
#' @import BART
#'
#' @export
simulate = function(N = 3000, discretize =c(3, 5, 10), setting = c(1,2,3), seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  Z = sapply(c(0.7, 0.5, 0.4), function(p){
    return(rbinom(N, 1, p))
  })
  X = rnorm(N)
  W = rnorm(N)
  auX = sapply(discretize, function(q) gtools::quantcut(X, q = q, labels = FALSE))
  colnames(auX) = paste0("auX_", discretize)
  auW = sapply(discretize, function(q) gtools::quantcut(W, q = q, labels = FALSE))
  colnames(auW) = paste0("auW_", discretize)

  if(setting == 1){
    pi = rstanarm::invlogit(- 1.25 - Z[,1] + 1.25*Z[,2] - 0.75*Z[,3] + 0.75*X - 0.10*X^2)
    Y1 = rnorm(N, 15 + 2.5*Z[,1] - Z[,2] + Z[,3] - 2*X + 3.75*X^2, 3)
    Y2 = rbinom(N, 1, rstanarm::invlogit(-2.5 + 0.75*Z[,1] - 2.5*Z[,2] + 1.5*Z[,3] - 0.25*X + 1.5*X^2))
  }
  if(setting == 2){
    pi = rstanarm::invlogit(- 1.25 - Z[,1] + 1.25*Z[,2] - 0.75*Z[,3] + 0.75*W - 0.10*W^2)
    Y1 = rnorm(N, 15 + 2.5*Z[,1] - Z[,2] + Z[,3] - 2*X + 3.75*X^2, 3)
    Y2 = rbinom(N, 1, rstanarm::invlogit(-2.5 + 0.75*Z[,1] - 2.5*Z[,2] + 1.5*Z[,3] - 0.25*X + 1.5*X^2))
  }
  if(setting == 3){
    pi = invlogit(-0.9 - 0.5 * Z[,1] + 0.75 * Z[,2] - Z[,3] + 0.5 * X - 0.1 * X^2 + 0.5 * Z[,1] * X + 0.25 * Z[,1] * X^2)
    Y1 = rnorm(N, 15 + 2.5 * Z[,1] - Z[,2] + Z[,3] - 2 * X + X^2 + Z[,1] * X - 2 * Z[,1] * X^2, 2)
    Y2 = rbinom(N, 1, invlogit(-1.75 + 0.75 * Z[,1] - 1.5 * Z[,2] + 1.5 * Z[,3] - 1.5 * X + X^2 + Z[,1] * X - 2 * Z[,1] * X^2))
  }
  if(setting == 4){
    pi = invlogit(-0.9 - 0.5 * Z[,1] + 0.75 * Z[,2] - Z[,3] + 0.5 * X - 0.1 * X^2 + 0.5 * Z[,1] * X + 0.75 * Z[,1] * X^2)
    Y1 = rnorm(N, 15 + 2.5 * Z[,1] - Z[,2] + Z[,3] - 2 * X + X^2 + 1 * Z[,1] * X - 4 * Z[,1] * X^2, 2)
    Y2 = rbinom(N, 1, invlogit(-1.75 + 0.75 * Z[,1] - 1.5 * Z[,2] + 1.5 * Z[,3] - 1.5 * X + X^2 + Z[,1] * X - 2.5 * Z[,1] * X^2))
  }
  population = data.frame(
    id = seq(1, N),
    Z = Z,
    X = X,
    W = W,
    auX = auX,
    auW = auW,
    Y1 = Y1,
    Y2 = Y2,
    true_pi = pi,
    logit_true_pi = rstanarm::logit(pi)
  )

  U = runif(N)
  population$inclusion = ifelse(pi > U, TRUE, FALSE)
  samples = population[population$inclusion == TRUE,]
  colnames(population) = c("id", "Z1", "Z2", "Z3", "X", "W", colnames(auX), colnames(auW), "Y1", "Y2", "true_pi", "logit_true_pi", "inclusion")
  colnames(samples) = c("id", "Z1", "Z2", "Z3", "X", "W", colnames(auX), colnames(auW), "Y1", "Y2", "true_pi", "logit_true_pi", "inclusion")

  ps_model = BART::pbart(as.matrix(population[, c("Z1", "Z2", "Z3", "X")]), population$inclusion, ntree=50, keepevery= 1, ndpost = 100)
  population$estimated_pi = predict(ps_model, as.matrix(population[, c("Z1", "Z2", "Z3", "X")]))$prob.test.mean
  samples$estimated_pi = predict(ps_model, as.matrix(samples[, c("Z1", "Z2", "Z3", "X")]))$prob.test.mean
  population$logit_estimated_pi = logit(population$estimated_pi)
  samples$logit_estimated_pi = logit(samples$estimated_pi)

  ps_model_with_W <- BART::pbart(as.matrix(population[, c("Z1", "Z2", "Z3", "X", "W")]), population$inclusion, ntree=50, keepevery= 1, ndpost = 100)
  population$estimated_pi_with_W = predict(ps_model_with_W, as.matrix(population[, c("Z1", "Z2", "Z3", "X", "W")]))$prob.test.mean
  samples$estimated_pi_with_W = predict(ps_model_with_W, as.matrix(samples[, c("Z1", "Z2", "Z3", "X", "W")]))$prob.test.mean
  population$logit_estimated_pi_with_W = logit(population$estimated_pi_with_W)
  samples$logit_estimated_pi_with_W = logit(samples$estimated_pi_with_W)

  samples$inclusion = NULL
  population$inclusion = NULL
  return(list(population = dplyr::as_tibble(population), samples = dplyr::as_tibble(samples)))
}



#' Weighted or Unweighted Sample Mean
#'
#' This function estimates the sample mean of an outcome variable using either weighted or unweighted
#' methods. It supports calculating the sample mean with finite population correction (FPC) when a
#' population dataset is provided. The method can also compute confidence intervals (CIs) for the sample mean
#' using the specified distribution family (Gaussian or Binomial).
#'
#' @param svysmpl A dataframe or tibble representing the sample data (`samples`). This should contain
#'                the outcome variable and any additional covariates.
#' @param svyVar The outcome variable to estimate the sample mean for (e.g., `Y1`).
#' @param svypopu A dataframe or tibble representing the population data (`population`). This is used to compute
#'                the finite population correction (FPC) when calculating the sample mean. Default is `NULL`.
#' @param subset A character vector representing filtering conditions to select subsets of the sample and population.
#'               Default is `NULL`, in which case the analysis is performed on the entire dataset. If subsets are specified,
#'               estimates for both the whole data and the subsets will be calculated.
#' @param family The distribution family of the outcome variable. Supported options are:
#'               \code{\link[stats]{gaussian}} for continuous outcomes and \code{\link[stats]{binomial}} for binary outcomes.
#' @param invlvls A numeric vector specifying the confidence levels (CIs) for the estimators. If more than
#'                one value is provided, multiple CIs will be calculated.
#' @param weights A numeric vector of case weights. The length should match the number of cases in `svysmpl`.
#'                These weights are used for calculating the weighted sample mean.
#'
#' @return A list, where each element contains the sample mean estimate and CIs for a subset or the entire data. The list includes:
#'         - `est`: The sample mean estimate.
#'         - `se`: The standard error of the sample mean estimate.
#'         - `tCI`: The confidence intervals for the sample mean.
#'         - `sample_size`: The sample size for the subset or entire dataset.
#'         - `population_size`: The population size, if a population dataset is provided (applicable to finite population correction).
#'         The list is returned for each subset specified.
#'
#' @examples
#' ## Simulate data with nonlinear association (setting 3).
#' data = simulate(N = 3000, discretize = 3, setting = 3, seed = 123)
#' population = data$population  # Population data (3000 cases)
#' samples = data$samples        # Sample data (600 cases)
#' ipw = 1 / samples$true_pi    # Compute inverse probability weights
#'
#' ## Estimate the weighted sample mean with IPW
#' IPW_sample_mean = uwt(svysmpl = samples, svyVar = "Y1", svypopu = population,
#'                       subset = c("Z1 == 1 & Z2 == 1"), family = gaussian(),
#'                       invlvls = c(0.95), weights = ipw)
#' IPW_sample_mean
#'
#' ## Estimate the unweighted sample mean
#' unweighted_sample_mean = uwt(svysmpl = samples, svyVar = "Y1", svypopu = population,
#'                               subset = NULL, family = gaussian(), invlvls = c(0.95), weights = NULL)
#' unweighted_sample_mean
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import rlang
#'
#' @export

uwt <- function(svysmpl, svyVar, svypopu = NULL, subset = NULL, family = gaussian(), invlvls, weights = NULL) {
  subset = c("T", subset)
  if(is.null(svypopu)){
    if(is.null(weights))
      des <- survey::svydesign(ids = ~1, weights = ~1, data = svysmpl)
    else{
      des <- survey::svydesign(ids = ~1, weights = ~weights, data = svysmpl)
    }
  }else{
    message("population parameter is specified, so the finite population correction will be calculated for sample mean.\n")
    svysmpl$fpc <- nrow(svypopu)
    if(is.null(weights))
      des <- survey::svydesign(ids = ~1, weights = ~1, data = svysmpl, fpc = ~fpc)
    else{
      des <- survey::svydesign(ids = ~1, weights = ~weights, data = svysmpl, fpc = ~fpc)
    }
  }
  infr = sapply(subset, function(s){
    des = subset(des, eval(parse(text = s)))

    if(family$family == "binomial"){
      suppressWarnings(desc <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula() %>%
                                        survey::svyciprop(des, method = "logit", level = lv), simplify = FALSE))
      tCI <- lapply(desc, function(i){
        ci = confint(i, df = survey::degf(des), parm = svyVar)
        colnames(ci) = stringr::str_replace(colnames(ci), "%", " %")
        ci
      })

      tCI = do.call("cbind", tCI)
      desc = desc[[1]]
    }
    if(family$family == "gaussian"){
      desc <- paste0('~', svyVar) %>% as.formula() %>% survey::svymean(des)
      tCI <- sapply(invlvls, confint, object = desc, parm = svyVar, simplify = FALSE)
      tCI = do.call("cbind", tCI)
    }

    # get estmates and standard error
    if(!is.null(svypopu)){
      infr <- cbind(est = desc[svyVar], se = sqrt(diag(vcov(desc))), tCI, sample_size = survey::degf(des) + 1, population_size =  nrow(dplyr::filter(svypopu, eval(parse(text = s)))))
    }else{
      infr <- cbind(est = desc[svyVar], se = sqrt(diag(vcov(desc))), tCI, sample_size = survey::degf(des) + 1)
    }
    if(is.null(weights))
      rownames(infr) = "sample_mean"
    else{
      rownames(infr) = "Weighted-sample_mean"
    }
    infr
  }, simplify = FALSE)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}

#' Weighted or Unweighted Raking Estimator
#'
#' This function estimates the weighted or unweighted raking adjustment for survey data. Raking adjusts the sample
#' weights to match the marginal distributions of auxiliary variables in the population. It supports both weighted
#' and unweighted estimations for a variety of outcome variables, including Gaussian (continuous) and Binomial (binary) outcomes.
#'
#' @param svysmpl A dataframe or tibble representing the sample data (`samples`). This should contain the outcome variable
#'                and any auxiliary variables.
#' @param svypopu A dataframe or tibble representing the population data (`population`). This is used to compute the
#'                finite population correction (FPC) for raking.
#' @param auxVars A character vector containing the names of auxiliary variables to be used for raking. These variables will
#'                be used to adjust the weights.
#' @param svyVar The outcome variable for which the raking estimate is calculated.
#' @param subset A character vector representing filtering conditions to select subsets of the sample and population.
#'               Default is `NULL`, in which case the analysis is performed on the entire dataset. If subsets are specified,
#'               estimates for both the whole data and the subsets will be calculated.
#' @param family The distribution family of the outcome variable. Supported options are:
#'               \code{\link[stats]{gaussian}} for continuous outcomes and \code{\link[stats]{binomial}} for binary outcomes.
#' @param invlvls A numeric vector specifying the confidence levels for the raking estimators. If more than one value is
#'                provided, multiple CIs will be calculated.
#' @param weights A numeric vector of case weights. The length should match the number of cases in `svysmpl`.
#'                These weights are used in the weighted raking adjustment.
#' @param maxiter An integer specifying the maximum number of iterations for the raking algorithm. Default is 50.
#'
#' @return A list where each element contains the raking estimate and confidence intervals (CIs) for a subset or the entire
#'         dataset. The list includes:
#'         - `est`: The raking estimate for the outcome variable.
#'         - `se`: The standard error of the estimate.
#'         - `tCI`: Confidence intervals for the estimate.
#'         - `sample_size`: The sample size for the subset or entire dataset.
#'         - `population_size`: The population size, if provided, including the finite population correction (FPC).
#'
#' @examples
#' ## Simulate data with nonlinear association (setting 3).
#' data = simulate(N = 3000, discretize = 3, setting = 3, seed = 123)
#' population = data$population  # Population data (3000 cases)
#' samples = data$samples        # Sample data (600 cases)
#' ipw = 1 / samples$true_pi    # Compute inverse probability weights
#'
#' ## Perform weighted raking with auxiliary variables
#' auxVars = c("Z1", "Z2", "Z3")
#' Weighted_rake = rake_wt(svysmpl = samples, svypopu = population, auxVars = auxVars,
#'                         svyVar = "Y1", subset = NULL, family = gaussian(),
#'                         invlvls = c(0.95), weights = ipw, maxiter = 50)
#' Weighted_rake
#'
#' ## Perform unweighted raking
#' Unweighted_rake = rake_wt(svysmpl = samples, svypopu = population, auxVars = auxVars,
#'                           svyVar = "Y1", subset = NULL, family = gaussian(),
#'                           invlvls = c(0.95), weights = NULL, maxiter = 50)
#' Unweighted_rake
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
#' @export

rake_wt <- function(svysmpl, svypopu, auxVars, svyVar, subset = NULL, family = gaussian(), invlvls, weights = NULL, maxiter = 50) {
  subset = c("T", subset)
  svysmpl$fpc <- nrow(svypopu)
  if(is.null(weights))
    des <- survey::svydesign(ids = ~1, weights = ~1, data = svysmpl, fpc = ~fpc)
  else{
    des <- survey::svydesign(ids = ~1, weights = ~weights, data = svysmpl, fpc = ~fpc)
  }
  popu_tab <- lapply(auxVars, function (Var) {
    Var = paste(all.vars(as.formula(paste0("~", Var))), collapse = "+")
    fmla <- paste('~', Var, sep = ' ') %>% as.formula()
    tab <- xtabs(fmla, svypopu)
    return(tab)
  } )
  fmla_list <- lapply(auxVars, function (Var) {
    Var = paste(all.vars(as.formula(paste0("~", Var))), collapse = "+")
    paste('~', Var, sep = ' ') %>% as.formula() %>% return()
  } )
  rakingobj <- survey::rake(des, fmla_list, popu_tab,
                            control = list(maxit = maxiter, epsilon = 1, verbose = FALSE))
  infr = sapply(subset, function(s){
    rakingobj = subset(rakingobj, eval(parse(text = s)))
    if(family$family == "binomial"){
      suppressWarnings(rakest <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula() %>%
                                          survey::svyciprop(rakingobj, method = "logit", level = lv), simplify = FALSE))
      tCI <- lapply(rakest, function(i){
        ci = confint(i, df = survey::degf(rakingobj), parm = svyVar)
        colnames(ci) = stringr::str_replace(colnames(ci), "%", " %")
        ci
      })
      tCI = do.call("cbind", tCI)
      rakest = rakest[[1]]
    }
    if(family$family == "gaussian"){
      rakest <- paste0('~', svyVar) %>% as.formula() %>%
        survey::svymean(rakingobj)
      tCI <- sapply(invlvls, confint, object = rakest, df = survey::degf(rakingobj),
                    parm = svyVar, simplify = FALSE)
      tCI = do.call("cbind", tCI)
    }
    infr <- cbind(est = rakest[svyVar], se = sqrt(diag(vcov(rakest))), tCI, sample_size = survey::degf(rakingobj) + 1, population_size = nrow(dplyr::filter(svypopu, eval(parse(text = s)))))
    if(is.null(weights))
      rownames(infr) = "rake"
    else{
      rownames(infr) = "Weighted-rake"
    }
    return(infr)
  }, simplify = FALSE)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}

#' Weighted or Unweighted Post-Stratification Estimator
#'
#' This function performs post-stratification adjustment for survey data, which adjusts the sample
#' weights to match the marginal distributions of auxiliary variables in the population. It supports
#' both weighted and unweighted estimations for various outcome variables, including Gaussian (continuous)
#' and Binomial (binary) outcomes. The function computes estimates and confidence intervals (CIs) for
#' the outcome variable using post-stratification based on the specified auxiliary variables.
#'
#' @param svysmpl A dataframe or tibble representing the sample data (`samples`). This should contain
#'                the outcome variable and any auxiliary variables.
#' @param svypopu A dataframe or tibble representing the population data (`population`). This is used to compute the
#'                finite population correction (FPC) for post-stratification.
#' @param auxVars A character vector containing the names of auxiliary variables to be used for post-stratification.
#'                These variables will be used to adjust the weights.
#' @param svyVar The outcome variable for which the post-stratification estimate is calculated.
#' @param subset A character vector representing filtering conditions to select subsets of the sample and population.
#'               Default is `NULL`, in which case the analysis is performed on the entire dataset. If subsets are specified,
#'               estimates for both the whole data and the subsets will also be calculated.
#' @param family The distribution family of the outcome variable. Supported options are:
#'               \code{\link[stats]{gaussian}} for continuous outcomes and \code{\link[stats]{binomial}} for binary outcomes.
#' @param invlvls A numeric vector specifying the confidence levels for the post-stratification estimators.
#'                If more than one value is provided, multiple CIs will be calculated.
#' @param weights A numeric vector of case weights. The length should match the number of cases in `svysmpl`.
#'                These weights are used in the weighted post-stratification adjustment.
#'
#' @return A list where each element contains the post-stratification estimate and confidence intervals (CIs) for a subset or the entire
#'         dataset. The list includes:
#'         - `est`: The post-stratification estimate for the outcome variable.
#'         - `se`: The standard error of the estimate.
#'         - `tCI`: The confidence intervals for the estimate.
#'         - `sample_size`: The sample size for the subset or entire dataset.
#'         - `population_size`: The population size, if provided, including the finite population correction (FPC).
#'
#' @examples
#' ## Simulate data with nonlinear association (setting 3).
#' data = simulate(N = 3000, discretize = 3, setting = 3, seed = 123)
#' population = data$population  # Population data (3000 cases)
#' samples = data$samples        # Sample data (600 cases)
#' ipw = 1 / samples$true_pi    # Compute inverse probability weights
#'
#' ## Perform weighted post-stratification with auxiliary variables
#' auxVars = c("Z1", "Z2", "Z3")
#' Weighted_postStratify = postStr_wt(svysmpl = samples, svypopu = population, auxVars = auxVars,
#'                                    svyVar = "Y1", subset = NULL, family = gaussian(),
#'                                    invlvls = c(0.95), weights = ipw)
#' Weighted_postStratify
#'
#' ## Perform unweighted post-stratification
#' Unweighted_postStratify = postStr_wt(svysmpl = samples, svypopu = population, auxVars = auxVars,
#'                                      svyVar = "Y1", subset = NULL, family = gaussian(),
#'                                      invlvls = c(0.95), weights = NULL)
#' Unweighted_postStratify
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
#' @export

postStr_wt <- function(svysmpl, svypopu, auxVars, svyVar, subset = NULL, family = gaussian(), invlvls, weights = NULL) {
  subset = c("T", subset)
  svysmpl$fpc <- nrow(svypopu) # get population sample size

  # specify svydesign
  if(is.null(weights))
    des <- survey::svydesign(ids = ~1, weights = ~1, data = svysmpl, fpc = ~fpc)
  else{
    des <- survey::svydesign(ids = ~1, weights = ~weights, data = svysmpl, fpc = ~fpc)
  }


  fmla <- paste(auxVars, collapse = "+")
  fmla = paste('~', fmla, sep = ' ') %>% as.formula()

  # set for post-stratification
  tab <- xtabs(fmla, svypopu)
  PSobj <- survey::postStratify(des, fmla, tab, partial = TRUE)
  #d = population[population$Z3 == 0 &population$auX_5  == 1,  ]
  #table(d$Z1, d$Z2)
  infr = sapply(subset, function(s){

    # make auxVar as a readable formula, for example, if auxVar = c(Z1, Z2, Z3), the result formula is ~Z1 + Z2 + Z3
    PSobj = subset(PSobj, eval(parse(text = s)))
    #svytable(fmla, PSobj, round = TRUE)

    # get estimates and confidence intervals
    # this function allows users specify multiple confidence levels, such as invlvls = c(0.95, 0.8)
    # so sapply function will calculate every confidence intervals separately
    if(family$family == "binomial"){
      suppressWarnings(PSest <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula() %>%
                                         survey::svyciprop(PSobj, method = "logit", level = lv), simplify = FALSE))
      tCI <- lapply(PSest, function(i){
        ci = confint(i, df = survey::degf(PSobj), parm = svyVar)
        colnames(ci) = str_replace(colnames(ci), "%", " %")
        ci
      })
      tCI = do.call("cbind", tCI)
      PSest = PSest[[1]]
    }
    if(family$family == "gaussian"){
      PSest <- paste0('~', svyVar) %>% as.formula() %>% survey::svymean(PSobj)
      tCI <- sapply(invlvls, confint, object = PSest, df = survey::degf(PSobj),
                    parm = svyVar, simplify = FALSE)
      tCI = do.call("cbind", tCI)
    }
    # get estmates and standard error
    infr <- cbind(est = PSest[svyVar], se = sqrt(diag(vcov(PSest))), tCI, sample_size = survey::degf(PSobj) + 1, population_size = nrow(dplyr::filter(svypopu, eval(parse(text = s)))))
    if(is.null(weights))
      rownames(infr) = "postStratify"
    else{
      rownames(infr) = "Weighted-postStratify"
    }
    return(infr)
  }, simplify = FALSE)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}


#' Bayesian Survey Model Estimation
#'
#' This function fits a Bayesian model using Stan for survey data. It allows you to specify
#' the outcome formula, the function for Stan, and apply different types of survey analysis, including
#' weighted or unweighted models, for both sample and population data. The function supports
#' posterior estimation, confidence intervals (CIs), and MCMC diagnostics.
#'
#' @param svysmpl A dataframe or tibble representing the sample data (`samples`). This should contain
#'                the outcome variable and any additional covariates.
#' @param svypopu A dataframe or tibble representing the population data (`population`).
#'                This should contain all variables in the model.
#' @param outcome_formula A formula for Stan, specifying the outcome and predictors in the model.
#' @param BayesFun The name of the Stan function to be used for fitting the Bayesian model.
#' @param subset A character vector representing filtering conditions to select subsets of the sample and population.
#'               Default is `NULL`, in which case the analysis is performed on the entire dataset. If specified,
#'               estimates for both the whole data and the subsets will be calculated.
#' @param family The distribution family for the outcome variable. Currently, the following options are supported:
#'               \code{\link[stats]{gaussian}} for continuous outcomes and \code{\link[stats]{binomial}} for binary outcomes.
#' @param invlvls A numeric vector specifying the confidence levels for the credible intervals (CIs). If more than
#'                one value is specified, multiple CIs will be calculated.
#' @param weights A numeric vector of case weights. The length of this vector should match the number of cases in `svysmpl`.
#'                These weights will be used in the Bayesian model for weighted estimation.
#' @param nskip An integer specifying the number of burn-in iterations for each chain in the MCMC for Stan models.
#'              Default is 1000.
#' @param npost An integer specifying the number of posterior sampling iterations for each chain in the MCMC for Stan models.
#'              Default is 1000.
#' @param nchain An integer specifying the number of MCMC chains for Stan models. Default is 4.
#' @param printmod A logical scalar; if `TRUE`, posterior estimates will be printed.
#' @param doFigure A logical scalar; if `TRUE`, MCMC diagnostic plots will be generated.
#' @param useTrueSample A logical scalar; if `TRUE`, the estimator will use true sample information.
#' @param stan_verbose A logical scalar; if `TRUE`, MCMC information will be printed during Stan model fitting.
#' @param HPD_CI A logical scalar; if `TRUE`, the calculated credible intervals will be highest posterior density intervals (HPD).
#'                    Otherwise, symmetric intervals will be used. Default is `FALSE`.
#' @param seed An integer specifying the random seed for reproducibility. Default is `NULL`.
#'
#' @return A list containing the Bayesian estimates and confidence intervals (CIs) for each subset or the entire dataset.
#'         Each element in the list includes:
#'         - `estimate`: The Bayesian point estimate for the outcome.
#'         - `CI`: The credible intervals for the outcome estimate.
#'         - Other elements based on the specified confidence levels in `invlvls`.
#'
#' @examples
#' \donttest{
#' ## Example usage with survey data:
#' ## Simulate sample and population data
#' data = simulate(N = 3000, discretize = 3, setting = 3, seed = 123)
#' population = data$population  # Get population data
#' samples = data$samples        # Get sample data
#' ipw = 1 / samples$true_pi    # Compute inverse probability weights
#'
#' ## Define outcome formula and Stan function
#' outcome_formula = "Y1 ~ Z1 + Z2 + Z3 + (1|auX_3)"
#' BayesFun = "stan_glmer"
#'
#' ## Fit Bayesian model using weighted survey data
#' bayes_model = svyBayesmod(svysmpl = samples, svypopu = population,
#'                           outcome_formula = outcome_formula,
#'                           BayesFun = BayesFun, weights = ipw,
#'                           family = gaussian(), nskip = 2000, npost = 2000,
#'                           nchain = 2, printmod = TRUE, invlvls = 0.95, stan_verbose = TRUE)
#' }
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
#' @export

svyBayesmod <- function(svysmpl, svypopu, outcome_formula, BayesFun, subset = NULL, family = gaussian(), invlvls, weights = NULL, nskip = 1000, npost = 1000, nchain = 4, printmod = TRUE, doFigure = FALSE, useTrueSample = FALSE, stan_verbose = FALSE, HPD_CI = FALSE, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  #print(outcome_formula)
  subset = c("T", subset)
  fmla <- outcome_formula[1] %>% as.formula()
  #print("start fit model")
  if(!is.null(seed)){
    if(is.null(weights)){
      if(!is.na(outcome_formula[2])){
        bayesmod <- get(BayesFun)(fmla, family = family, random = as.formula(outcome_formula[2]), data = svysmpl, iter = nskip + npost, warmup = nskip, chains = nchain, refresh = stan_verbose, seed = seed)
      }else{
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, chains = nchain, refresh = stan_verbose, seed = seed)
      }
    }
    else{
      if(!is.na(outcome_formula[2])){
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, random = as.formula(outcome_formula[2]), refresh = stan_verbose, seed = seed)
      }else{
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, refresh = stan_verbose, seed = seed)
      }
      #bayesmod <- get(BayesFun)(fmla, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, random = as.formula(formula[2]))
    }
  }else{
    if(is.null(weights)){
      if(!is.na(outcome_formula[2])){
        bayesmod <- get(BayesFun)(fmla, family = family, random = as.formula(outcome_formula[2]), data = svysmpl, iter = nskip + npost, warmup = nskip, chains = nchain, refresh = stan_verbose)
      }else{
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, chains = nchain, refresh = stan_verbose)
      }
    }
    else{
      if(!is.na(outcome_formula[2])){
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, random = as.formula(outcome_formula[2]), refresh = stan_verbose)
      }else{
        bayesmod <- get(BayesFun)(fmla, family = family, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, refresh = stan_verbose)
      }
      #bayesmod <- get(BayesFun)(fmla, data = svysmpl, iter = nskip + npost, warmup = nskip, weights = weights, chains = nchain, random = as.formula(formula[2]))
    }
  }

  #print("end fit model")
  #if (meth %in% c('Cov', 'catCov') )  bayesmod <- stan_glm(fmla, data = svysmpl)
  #if (meth %in% c('GAMcatCov' , 'GAM-Cov', 'GAM-PS', 'GAM-lgtPS', 'GAM-slgtPS') ) bayesmod <- stan_gamm4(fmla, data = svysmpl)
  #bayesmod <- update(bayesmod, iter = 500)

  if (printmod) summary(bayesmod) %>% print()

  if (doFigure) {
    p1 <- pp_check(bayesmod)
    p2 <- pp_check(bayesmod, plotfun = 'scatter_avg')
    p3 <- pp_check(bayesmod, plotfun = 'stat_2d', stat = c('mean', 'sd'))
    p4 <- pp_check(bayesmod, plotfun = "error_binned")
    ppcFig <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
    ppcFig <- gridExtra::grid.arrange(ppcFig, p4, ncol = 1)
    ppcFig
  }
  if(useTrueSample == FALSE){
    infr = sapply(subset, function(s){
      #print(s)
      svypopu1 = dplyr::filter(svypopu, eval(parse(text = s)))
      yhats <- posterior_epred(bayesmod, svypopu1)

      post_est = rowMeans(yhats)

      tCI = sapply(invlvls, function(level){
        if(HPD_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, probs = level, names = TRUE)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = TRUE)
          names(ci) = str_replace(names(ci), "%", " %")
        }

        ci
      }, simplify = FALSE)
      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(dplyr::filter(svysmpl, eval(parse(text = s)))), population_size = nrow(svypopu1)))
    }, simplify = FALSE)
    names(infr)[1] = "All"
    return(infr)
  }else{
    infr = sapply(subset, function(s){
      svypopu1 = dplyr::filter(svypopu, eval(parse(text = s)))
      svysmpl1 = dplyr::filter(svysmpl, eval(parse(text = s)))

      yhats_pop <- posterior_predict(bayesmod, svypopu1, re.form = NA, seed = seed) # npost, non-sample_size
      yhats_sam <- posterior_predict(bayesmod, svysmpl1, re.form = NA, seed = seed)

      yhats_pop_tot <- yhats_pop %>% apply(1, sum) # npost * 1
      yhats_sam_tot <- yhats_sam %>% apply(1, sum) # npost * 1
      yobs_tot <- sum(svysmpl1 %>% dplyr::select(stringr::str_trim(stringr::str_split_1(outcome_formula[1], "~")[1]))) # a single value
      post_est = (yhats_pop_tot + yobs_tot - yhats_sam_tot) / nrow(svypopu1) #npost
      tCI = sapply(invlvls, function(level){
        #confint(post_est, level = 0.95)
        if(HPD_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, prob = level, names = TRUE)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = TRUE)
          names(ci) = str_replace(names(ci), "%", " %")
        }
        ci
      }, simplify = FALSE)

      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(svysmpl1), population_size = nrow(svypopu1)))
      return(infr)
    }, simplify = FALSE)
    names(infr)[1] = "All"
    return(infr)
  }
}



#' Auxiliary Variables in Survey Analysis
#'
#' This function provides a user-friendly interface for various estimators in survey analysis
#' when working with discretized auxiliary variables. Probability surveys often use continuous
#' data from administrative records as auxiliary variables, but the utility of this data is
#' diminished when discretized for confidentiality purposes. This package offers different estimators
#' that handle discretized auxiliary variables effectively.
#'
#' The available estimators include:
#' - Weighted or unweighted sample mean
#' - Weighted or unweighted raking
#' - Weighted or unweighted post-stratification
#' - Bayesian methods:
#'   - BART (Bayesian Additive Regression Trees)
#'   - MRP (Multilevel Regression with Poststratification)
#'   - GAMP (Generalized Additive Model of Response Propensity)
#'   - Weighted linear regression
#'
#' These Bayesian models are implemented using the **rstan** and **rstanarm** packages.
#'
#' @param formula A string or formula specifying the outcome model. For non-model-based methods
#' (e.g., sample mean, raking, post-stratification), only include the outcome variable (e.g., "~Y").
#' For model-based methods (e.g., MRP, GAMP, linear regression), additional fixed effect predictors can
#' be specified, such as "Y ~ X1 + X2 + I(X^2)". For GAMP, smooth functions can be specified as
#' "Y ~ X1 + s(X2, 10) + s(X3, by = X1)". Categorical variables are automatically treated as dummy variables
#' in model-based methods.
#' @param auxiliary A string specifying the formula for the auxiliary variables. For sample mean and
#' BART, this should be `NULL`. For raking, post-stratification, and GAMP, this should be an additive model
#' (e.g., "Z1 + Z2 + Z3"). For MRP, specify random effects for terms in this parameter, such as "Z1 + Z2 + Z3"
#' or "Z1 + Z2:Z3".
#' @param samples A dataframe or tibble containing all variables specified in `formula` and `auxiliary`.
#' This is typically a subset of the population.
#' @param population A dataframe or tibble containing all variables specified in `formula` and `auxiliary`.
#' This is the entire population used for estimation.
#' @param subset A character vector representing filtering conditions to select subsets of `samples` and `population`.
#' Default is `NULL`, in which case the analysis is performed on the entire dataset. If subsets are specified,
#' estimates for both the whole data and the subsets will be calculated.
#' @param family The distribution family of the outcome variable. Supported options are:
#' \code{\link[stats]{gaussian}} for continuous outcomes and \code{\link[stats]{binomial}} for binary outcomes.
#' @param method A string specifying the model to use. Options include "sample_mean", "rake", "postStratify",
#' "MRP", "GAMP", "linear", and "BART".
#' @param weights A numeric vector of case weights. The length should match the number of cases in `samples`.
#' @param levels A numeric vector specifying the confidence levels for the confidence intervals (CIs).
#' Multiple values can be specified to calculate multiple CIs.
#' @param stan_verbose A logical scalar; if `TRUE`, prints all messages when running Stan models. Default is `FALSE`.
#' This parameter only applies to Bayesian models.
#' @param show_plot A logical scalar; if `TRUE`, shows diagnostic plots for Stan models. Default is `FALSE`.
#' This parameter only applies to Bayesian models.
#' @param nskip An integer specifying the number of burn-in iterations for each chain in MCMC for Stan models.
#' Default is `1000`. This parameter only applies to Bayesian models.
#' @param npost An integer specifying the number of posterior sampling iterations for each chain in MCMC for Stan models.
#' Default is `1000`. This parameter only applies to Bayesian models.
#' @param nchain An integer specifying the number of MCMC chains for Stan models. Default is `4`. This parameter
#' only applies to Bayesian models.
#' @param HPD_interval A logical scalar; if `TRUE`, calculates the highest posterior density (HPD) intervals for the
#' CIs of Stan models. Default is `FALSE`, in which case symmetric intervals are calculated. This parameter only applies
#' to Bayesian models.
#' @param seed An integer specifying the random seed for reproducibility. Default is `NULL`.
#'
#' @return A list containing the sample mean estimates and CIs for the subset and/or the whole dataset.
#'         Each element in the list includes:
#'         - `estimate`: The point estimate of the sample mean.
#'         - `CI`: Confidence intervals for the sample mean.
#'         - Other elements for each confidence level specified in `levels`.
#'
#' @examples
#' \donttest{
#' ## Simulate data with nonlinear association (setting 3).
#' data = simulate(N = 3000, discretize = 10, setting = 3, seed = 123)
#' population = data$population
#' samples = data$samples
#' ipw = 1 / samples$true_pi
#' true_mean = mean(population$Y1)
#'
#' ## IPW Sample Mean
#' IPW_sample_mean = auxsurvey("~Y1", auxiliary = NULL, weights = ipw,
#'                             samples = samples, population = population,
#'                             subset = c("Z1 == 1 & Z2 == 1"), method = "sample_mean",
#'                             levels = 0.95)
#'
#' ## Raking
#' rake = auxsurvey("~Y1", auxiliary = "Z1 + Z2 + Z3 + auX_10", samples = samples,
#'                  population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"),
#'                  method = "rake", levels = 0.95)
#'
#' ## MRP
#' MRP = auxsurvey("Y1 ~ 1 + Z1", auxiliary = "Z2 + Z3:auX_10", samples = samples,
#'                 population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"),
#'                 method = "MRP", levels = 0.95, nskip = 4000, npost = 4000,
#'                 nchain = 1, stan_verbose = FALSE, HPD_interval = TRUE)
#'
#' ## GAMP
#' GAMP = auxsurvey("Y1 ~ 1 + Z1 + Z2 + Z3", auxiliary = "s(auX_10) + s(logit_true_pi, by = Z1)",
#'                  samples = samples, population = population, method = "GAMP",
#'                  levels = 0.95, nskip = 4000, npost = 4000, nchain = 1,
#'                  stan_verbose = FALSE, HPD_interval = TRUE)
#'
#' ## BART
#' BART = auxsurvey("Y1 ~ Z1 + Z2 + Z3 + auX_10", auxiliary = NULL, samples = samples,
#'                  population = population, method = "BART", levels = 0.95,
#'                  nskip = 4000, npost = 4000, nchain = 1, HPD_interval = TRUE)
#'
#' }
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr mutate_at filter %>% pull
#' @import stringr
#' @import BART
#'
#' @export
auxsurvey <- function(formula, auxiliary = NULL, samples, population = NULL, subset = NULL, family = gaussian(), method = c("sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear", "BART"), weights = NULL, levels = c(0.95, 0.8, 0.5), stan_verbose = TRUE, show_plot = TRUE, nskip = 1000, npost = 1000, nchain = 4, HPD_interval = FALSE, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)

  svyVar = stringr::str_trim(stringr::str_split_1(as.character(formula), "~"))
  svyVar = svyVar[svyVar != ""][1]
  if(!is.null(auxiliary)){
    auxiliary = as.character(auxiliary)
    auxiliary = stringr::str_remove_all(auxiliary, " ")
    auxiliary = paste0(auxiliary, collapse = "+")
    if(!stringr::str_detect(auxiliary, "^~")){
      auxiliary = paste0("~",auxiliary)
    }
  }
  if(!is.null(population)){
    if(length(setdiff(union(all.vars(as.formula(auxiliary)), all.vars(as.formula(formula))), union(svyVar, colnames(population)))) > 0){
      stop(paste0("unidentified variables: ", paste0(setdiff(union(all.vars(as.formula(auxiliary)), all.vars(as.formula(formula))), colnames(population)), collapse = ", "), collapse = ", "))
    }
  }
  if(length(setdiff(union(all.vars(as.formula(auxiliary)), all.vars(as.formula(formula))), colnames(samples))) > 0){
    stop(paste0("unidentified variables: ", paste0(setdiff(union(all.vars(as.formula(auxiliary)), all.vars(as.formula(formula))), colnames(samples)), collapse = ", "), collapse = ", "))
  }

  method = match.arg(method)
  # if(is.null(propensity_score))
  #   weights = NULL
  # else{
  #   weights = 1 / propensity_score
  # }
  if(method == "sample_mean"){
    return(uwt(samples, svyVar, population, subset, family, levels, weights))
  }
  if(method == "rake"){
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = TRUE))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    auxiliary = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE))
    auxiliary = auxiliary[!is.na(auxiliary) & auxiliary != ""]
    return(rake_wt(samples, population, auxiliary, svyVar, subset, family = family, levels, weights, maxiter = 50))
  }
  if(method == "postStratify"){
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = TRUE))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    auxiliary = unique(c(all.vars(as.formula(auxiliary)), covariates))
    auxiliary = auxiliary[!is.na(auxiliary) & auxiliary != ""]
    return(postStr_wt(samples, population, auxiliary, svyVar, subset, family = family, levels, weights))
  }
  if(method == "MRP"){
    if(is.null(nskip)) nskip = 1000
    if(is.null(npost)) npost = 1000
    if(is.null(nchain)) nchain = 4
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = TRUE))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    samples = dplyr::mutate_at(samples, all.vars(as.formula(auxiliary)), as.factor)
    population = dplyr::mutate_at(population, all.vars(as.formula(auxiliary)), as.factor)
    auxiliary = stringr::str_replace_all(auxiliary, "\\*", ":")
    auxiliary = stringr::str_split_i(as.character(auxiliary), "~", 2)
    auxiliary = stringr::str_split(auxiliary, "\\+", simplify = TRUE)

    if(length(covariates) == 0){
      outcome_formula = paste0(paste0(svyVar, "~"), paste0("(1|", auxiliary, ")", collapse = "+"))
    }else{
      outcome_formula = paste(paste0(svyVar, "~", paste0(covariates, collapse = "+")), paste0("(1|", auxiliary, ")", collapse = "+"), sep= "+")
    }
    message("The formula for the MRP model is ", outcome_formula, "\n")
    MRP_est = svyBayesmod(samples, population, outcome_formula, "stan_glmer", subset, family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = show_plot, useTrueSample = FALSE, stan_verbose = stan_verbose, HPD_CI = HPD_interval)
    #return(MRP_est)
    MRP_est =  lapply(MRP_est, function(est){
      if(is.null(weights))
        rownames(est) = "MRP"
      else{
        rownames(est) = "Weighted-MRP"
      }
      return(est)
    })
    if(length(MRP_est) == 1)
      return(MRP_est[[1]])
    return(MRP_est)

  }
  if(method == "GAMP"){
    if(is.null(nskip)) nskip = 1000
    if(is.null(npost)) npost = 1000
    if(is.null(nchain)) nchain = 4
    #covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = TRUE))
    #if("." %in% covariates){
    #  covariates = setdiff(names(samples), svyVar)
    #}

    #samples_matrix = model.matrix(as.formula(paste(paste0("~", str_split_i(formula, "~", 2)), paste0(auxiliary, collapse = "+"), sep = "+")), data = samples)
    samples = dplyr::mutate_at(samples, all.vars(as.formula(auxiliary)), as.numeric)
    population = dplyr::mutate_at(population, intersect(all.vars(as.formula(auxiliary)), colnames(population)), as.numeric)

    auxiliary = stringr::str_replace_all(auxiliary, "\\*", ":")
    auxiliary_fixed = setdiff(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), all.vars(as.formula(auxiliary)))[stringr::str_detect(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), "s\\(.*\\)")]


    outcome_formula = paste(formula, paste(auxiliary_fixed, collapse = "+"), sep = "+")
    auxiliary_random = union(intersect(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), all.vars(as.formula(auxiliary))), stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE)[is.na(stringr::str_match(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), "s\\(.*\\)"))])

    if(length(auxiliary_random) > 0){
      samples = dplyr::mutate_at(samples, all.vars(as.formula(paste0("~", auxiliary_random))), as.factor)
      population = dplyr::mutate_at(population, all.vars(as.formula(paste0("~", auxiliary_random))), as.factor)
    }
    if(length(auxiliary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxiliary_random, ")", collapse = "+")))
      message("The formula for the GAMP model is ", paste0(outcome_formula[1], str_replace(outcome_formula[2], "~", "+"), collapse = ""), "\n")
    }else{
      outcome_formula = c(outcome_formula, NULL)
      message("The formula for the GAMP model is ", outcome_formula[1], "\n")
    }

    #outcome_formula = paste(formula, paste0(str_remove_all(auxiliary, "~"), collapse = "+"), paste0("(1|", auxiliary_random, ")", collapse = "+"),  sep = "+")

    GAMP_est = svyBayesmod(samples, population, outcome_formula, "stan_gamm4", subset, family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = FALSE, useTrueSample = TRUE, stan_verbose = stan_verbose, HPD_CI = HPD_interval, seed = seed)
    #print(GAMP_est)

    GAMP_est =  lapply(GAMP_est, function(est){
      if(is.null(weights))
        rownames(est) = "GAMP"
      else{
        rownames(est) = "Weighted-GAMP"
      }
      return(est)
    })
    if(length(GAMP_est) == 1)
      return(GAMP_est[[1]])
    return(GAMP_est)
  }
  if(method == "linear"){
    if(is.null(nskip)) nskip = 1000
    if(is.null(npost)) npost = 1000
    if(is.null(nchain)) nchain = 4
    samples = mutate_at(samples, all.vars(as.formula(auxiliary)), as.numeric)
    population = dplyr::mutate_at(population, intersect(all.vars(as.formula(auxiliary)), colnames(population)), as.numeric)


    auxiliary_fixed = setdiff(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), all.vars(as.formula(auxiliary)))


    outcome_formula = paste(formula, paste(auxiliary_fixed, collapse = "+"), collapse = "+")
    auxiliary_random = intersect(str_split(str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = TRUE), all.vars(as.formula(auxiliary)))
    samples = dplyr::mutate_at(samples, auxiliary_random, as.factor)
    population = dplyr::mutate_at(population, auxiliary_random, as.factor)
    if(length(auxiliary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxiliary_random, ")", collapse = "+")))
      message("The formula for the linear model is ", paste0(outcome_formula[1], str_replace(outcome_formula[2], "~", "+"), collapse = ""), "\n")
      Linear_est = svyBayesmod(samples, population, outcome_formula, "stan_glmer", family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = FALSE, useTrueSample = TRUE, stan_verbose = stan_verbose, HPD_CI = HPD_interval)
    }else{
      outcome_formula = c(outcome_formula, NULL)
      message("The formula for the linear model is ", outcome_formula, "\n")
      Linear_est = svyBayesmod(samples, population, outcome_formula, "stan_glm", subset, family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = FALSE, useTrueSample = TRUE, stan_verbose = stan_verbose, HPD_CI = HPD_interval)
    }


    Linear_est = lapply(Linear_est, function(est){
      if(is.null(weights))
        rownames(est) = "linear"
      else{
        rownames(est) = "Weighted-linear"
      }
      return(est)
    })
    if(length(Linear_est) == 1)
      return(Linear_est[[1]])

    #outcome_formula = paste(formula, paste0(str_remove_all(auxiliary, "~"), collapse = "+"), paste0("(1|", auxiliary_random, ")", collapse = "+"),  sep = "+")
    return(Linear_est)
  }
  if(method == "BART"){
    if(is.null(nskip)) nskip = 1000
    if(is.null(npost)) npost = 1000
    if(is.null(nchain)) nchain = 1
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = TRUE))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    for(i in covariates){
      if(is.character(samples[,i])){
        samples[,i] = as.factor(samples[,i])
        population[,i] = as.factor(population[,i])
      }
    }
    if(family$family == "binomial"){
      X_train = stats::model.matrix(as.formula(paste0("~", stringr::str_split_i(as.character(formula), "~", 2))), samples[, covariates])
      y_train = dplyr::pull(samples, svyVar)
      model <- BART::pbart(X_train, y_train, ndpost = npost, nskip = nskip, rm.const=FALSE)
    }
    if(family$family == "gaussian"){
      X_train = stats::model.matrix(as.formula(paste0("~", stringr::str_split_i(as.character(formula), "~", 2))), samples[, covariates])
      #X_train = samples[, covariates]
      y_train = dplyr::pull(samples, svyVar)
      model <- BART::wbart(X_train, y_train, ndpost = npost, nskip = nskip, rm.const=FALSE)
    }
    subset = c("T", subset)
    infr = sapply(subset, function(s){
      svypopu1 = dplyr::filter(population, eval(parse(text = s)))[, covariates]
      svysmpl1 = dplyr::filter(samples, eval(parse(text = s)))[, covariates]
      svypopu1 = stats::model.matrix(as.formula(paste0("~", stringr::str_split_i(as.character(formula), "~", 2))), svypopu1)
      svysmpl1 = stats::model.matrix(as.formula(paste0("~", stringr::str_split_i(as.character(formula), "~", 2))), svysmpl1)
      yhats_pop <- predict(model, svypopu1) # npost, non-sample_size
      yhats_sam <- predict(model, svysmpl1)
      if(family$family == "binomial"){
        yhats_pop = apply(yhats_pop$prob.test, c(1, 2), function(x) rbinom(1, size = 1, prob = x))
        yhats_sam = apply(yhats_sam$prob.test, c(1, 2), function(x) rbinom(1, size = 1, prob = x))
      }


      yhats_pop_tot <- yhats_pop %>% apply(1, sum) # npost * 1
      yhats_sam_tot <- yhats_sam %>% apply(1, sum) # npost * 1
      yobs_tot <- sum(dplyr::pull(dplyr::filter(samples, eval(parse(text = s))), svyVar)) # a single value
      post_est = (yhats_pop_tot + yobs_tot - yhats_sam_tot) / nrow(svypopu1) #npost
      tCI = sapply(levels, function(level){
        #confint(post_est, level = 0.95)
        if(HPD_interval){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, prob = level, names = TRUE)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = TRUE)
          names(ci) = stringr::str_replace(names(ci), "%", " %")
        }
        ci
      }, simplify = FALSE)

      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(svysmpl1), population_size = nrow(svypopu1)))
      return(infr)
    }, simplify = FALSE)
    names(infr)[1] = "All"
    BART_est =  lapply(infr, function(est){
      rownames(est) = "BART"
      return(est)
    })
    if(length(BART_est) == 1)
      return(BART_est[[1]])
    return(BART_est)
  }
}


