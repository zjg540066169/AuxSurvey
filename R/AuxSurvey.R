library(mgcv)
library(rstanarm)

#' Simulate
#'
#' @param N Number of population size
#' @param discretize Number of discretized categories
#' @param setting Simulation setting
#' @param seed Random seed
#'
#' @return simulated dataset
#' @export
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom gtools quantcut
#' @import mgcv
#' @importFrom dplyr filter %>% as_tibble
#' @import stringr
#' @import BART
#'
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
    pi = rstanarm::invlogit(-0.9 - 0.5 * Z[,1] + 0.75 * Z[,2] - Z[,3] + 0.5 * X - 0.05 * X^2 + 0.5 * Z[,1] * X - 0.75 * Z[,1] * X^2)
    Y1 = rnorm(N, 15 + 2.5 * Z[,1] - Z[,2] + Z[,3] - 2 * X + X^2 + Z[,1] * X - 2.5 * Z[,1] * X^2, 2)
    Y2 = rbinom(N, 1, rstanarm::invlogit(-1.75 + 0.75 * Z[,1] - 1.5 * Z[,2] + 1.5 * Z[,3] - 1.5 * X + X^2 + Z[,1] * X - 2.5 * Z[,1] * X^2))
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
  population$inclusion = ifelse(pi > U, T, F)
  samples = population[population$inclusion == T,]
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
#' @param svysmpl A dataframe or tibble to represent 'samples'.
#' @param svyVar Outcome variable.
#' @param svypopu population.
#' @param subset A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.
#' @param family The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.
#' @param invlvls A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.
#' @param weights A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.
#'
#' @return A list. Each element contains the sample mean estimate and CIs for a subset or the whole data analysis.
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import rlang
#'
uwt <- function(svysmpl, svyVar, svypopu = NULL, subset = NULL, family = gaussian(), invlvls, weights = NULL) {
  subset = c("T", subset)
  if(is.null(svypopu)){
    if(is.null(weights))
      des <- survey::svydesign(ids = ~1, weights = ~1, data = svysmpl)
    else{
      des <- survey::svydesign(ids = ~1, weights = ~weights, data = svysmpl)
    }
  }else{
    cat("population parameter is specified, so the finite population correction will be calculated for sample mean.\n")
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
                                        survey::svyciprop(des, method = "logit", level = lv), simplify = F))
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
      tCI <- sapply(invlvls, confint, object = desc, parm = svyVar, simplify = F)
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
  }, simplify = F)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}

#' Weighted or Unweighted Raking
#'
#' @param svysmpl A dataframe or tibble to represent 'samples'.
#' @param svypopu A dataframe or tibble to represent 'population'.
#' @param auxVars A character vector contains the names of auxiliary variables for raking.
#' @param svyVar Outcome variable.
#' @param subset A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.
#' @param family The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.
#' @param invlvls A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.
#' @param weights A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.
#' @param maxiter A integer to specify maximum iteration of raking.
#'
#' @return A list. Each element contains the raking estimate and CIs for a subset or the whole data analysis.
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
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
                                          survey::svyciprop(rakingobj, method = "logit", level = lv), simplify = F))
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
                    parm = svyVar, simplify = F)
      tCI = do.call("cbind", tCI)
    }
    infr <- cbind(est = rakest[svyVar], se = sqrt(diag(vcov(rakest))), tCI, sample_size = survey::degf(rakingobj) + 1, population_size = nrow(dplyr::filter(svypopu, eval(parse(text = s)))))
    if(is.null(weights))
      rownames(infr) = "rake"
    else{
      rownames(infr) = "Weighted-rake"
    }
    return(infr)
  }, simplify = F)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}

#' Weighted or Unweighted Post-Stratification
#'
#' @param svysmpl A dataframe or tibble to represent 'samples'.
#' @param svypopu A dataframe or tibble to represent 'population'.
#' @param auxVars A character vector contains the names of auxiliary variables for Post-Stratification.
#' @param svyVar Outcome variable.
#' @param subset A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.
#' @param family The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.
#' @param invlvls A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.
#' @param weights A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.
#'
#' @return A list. Each element contains the Post-Stratification estimate and CIs for a subset or the whole data analysis.
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
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
    #svytable(fmla, PSobj, round = T)

    # get estimates and confidence intervals
    # this function allows users specify multiple confidence levels, such as invlvls = c(0.95, 0.8)
    # so sapply function will calculate every confidence intervals separately
    if(family$family == "binomial"){
      suppressWarnings(PSest <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula() %>%
                                         survey::svyciprop(PSobj, method = "logit", level = lv), simplify = F))
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
                    parm = svyVar, simplify = F)
      tCI = do.call("cbind", tCI)
    }
    print(s)
    # get estmates and standard error
    infr <- cbind(est = PSest[svyVar], se = sqrt(diag(vcov(PSest))), tCI, sample_size = survey::degf(PSobj) + 1, population_size = nrow(dplyr::filter(svypopu, eval(parse(text = s)))))
    print(infr)
    if(is.null(weights))
      rownames(infr) = "postStratify"
    else{
      rownames(infr) = "Weighted-postStratify"
    }
    return(infr)
  }, simplify = F)
  names(infr)[1] = "All"
  if(length(infr) == 1)
    return(infr[[1]])
  return(infr)
}


#' Title
#'
#' @param svysmpl A dataframe or tibble to represent 'samples'.
#' @param svypopu A dataframe or tibble to represent 'population'.
#' @param outcome_formula A formula for stan.
#' @param BayesFun Function name for stan
#' @param subset A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.
#' @param family The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.
#' @param invlvls A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.
#' @param weights A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.
#' @param nskip A integer to specify the number of burn-in iterations of each chain in MCMC for stan models. Default is 1000.
#' @param npost A integer to specify the number of posterior sampling iteration of each chain in MCMC for stan models. Default is 1000.
#' @param nchain A integer to specify the number of MCMC chains for stan models. Default is 4.
#' @param printmod A logical to indicate if print posterior estimates
#' @param doFigure A logical to indicate if print MCMC figures
#' @param useTrueSample A logical to indicate if the estimator uses samples information
#' @param stan_verbose A logical to indicate if print MCMC information in stan
#' @param shortest_CI  A logical scalar; if true, the calculated credible intervals for stan models are highest posterior density intervals. Otherwise the intervals are symmetric. Default is false.
#' @return A list. Each element contains the Bayesian estimate and CIs for a subset or the whole data analysis.
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr filter %>%
#' @import stringr
#'
svyBayesmod <- function(svysmpl, svypopu, outcome_formula, BayesFun, subset = NULL, family = gaussian(), invlvls, weights = NULL, nskip = 1000, npost = 1000, nchain = 4, printmod = TRUE, doFigure = FALSE, useTrueSample = F, stan_verbose = F, shortest_CI = F) {
  #print(outcome_formula)
  subset = c("T", subset)
  fmla <- outcome_formula[1] %>% as.formula()
  #print("start fit model")
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
  if(useTrueSample == F){
    infr = sapply(subset, function(s){
      #print(s)
      svypopu1 = dplyr::filter(svypopu, eval(parse(text = s)))
      yhats <- posterior_epred(bayesmod, svypopu1)

      post_est = rowMeans(yhats)

      tCI = sapply(invlvls, function(level){
        if(shortest_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, probs = level, names = T)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = T)
          names(ci) = str_replace(names(ci), "%", " %")
        }

        ci
      }, simplify = F)
      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(dplyr::filter(svysmpl, eval(parse(text = s)))), population_size = nrow(svypopu1)))
    }, simplify = F)
    names(infr)[1] = "All"
    return(infr)
  }else{
    infr = sapply(subset, function(s){
      svypopu1 = dplyr::filter(svypopu, eval(parse(text = s)))
      svysmpl1 = dplyr::filter(svysmpl, eval(parse(text = s)))

      yhats_pop <- posterior_predict(bayesmod, svypopu1, re.form = NA) # npost, non-sample_size
      yhats_sam <- posterior_predict(bayesmod, svysmpl1, re.form = NA)

      yhats_pop_tot <- yhats_pop %>% apply(1, sum) # npost * 1
      yhats_sam_tot <- yhats_sam %>% apply(1, sum) # npost * 1
      yobs_tot <- sum(svysmpl1 %>% dplyr::select(stringr::str_trim(stringr::str_split_1(outcome_formula[1], "~")[1]))) # a single value
      post_est = (yhats_pop_tot + yobs_tot - yhats_sam_tot) / nrow(svypopu1) #npost
      tCI = sapply(invlvls, function(level){
        #confint(post_est, level = 0.95)
        if(shortest_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, prob = level, names = T)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = T)
          names(ci) = str_replace(names(ci), "%", " %")
        }
        ci
      }, simplify = F)

      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(svysmpl1), population_size = nrow(svypopu1)))
      return(infr)
    }, simplify = F)
    names(infr)[1] = "All"
    return(infr)
  }
}



#' Auxiliary Variables in Survey Analysis
#'@description Probability surveys often use auxiliary continuous data from administrative records, but the utility of this data is diminished when it is discretized for confidentiality.
#'
#'This function provides a user-friendly interface for different estimators with the discretized auxiliary variables.
#'
#'The estimators includes (weighted) sample mean, (weighted) raking, (weighted) post-stratification, and four Bayesian methods: BART, MRP, GAMP(Generalized additive model of response propensity) and (weighted) linear regression. The three Bayesian models are based on 'rstan' and 'rstanarm'.
#'
#' @param formula A string or formula for the specified formula for the outcome model. For non-model based methods(sample mean, raking, post-stratification), just include the outcome variable, such as "~Y". For model based methods (MRP, GAMP, LR), additional predictors can be specified as fixed effects term, such as "Y~X1+X2 + I(X^2)". For GAMP, smooth functions can be specified, such as "Y~X1 + s(X2, 10) + s(X3, by = X1)". Categorical variables are coded as dummy variables in model based methods.
#' @param auxiliary A string for the specified formula for the auxiliary variables. For sample mean and BART, just leave it as NULL. For raking, post-stratification, GAMP, use string for an additive model, such as "Z1 + Z2 + Z3". MRP specifies random effects for terms in this parameter, such as "Z1 + Z2 + Z3" or "Z1 + Z2:Z3".
#' @param samples A dataframe or tibble contains all variables in 'formula' and 'auxiliary'. This dataframe is a subset of 'population'.
#' @param population A dataframe or tibble contains all variables in 'formula' and 'auxiliary'.
#' @param subset A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.
#' @param family The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.
#' @param method A string specifying which model to use.
#' @param weights A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.
#' @param levels A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.
#' @param stan_verbose A logical scalar; if true, print all messages when running stan models. Default is false. This parameter only works for Bayesian models.
#' @param show_plot A logical scalar; if true, show some diagnostic plots for stan models. Default is false. This parameter only works for Bayesian models.
#' @param nskip An integer to specify the number of burn-in iterations of each chain in MCMC for stan models. Default is 1000. This parameter only works for Bayesian models.
#' @param npost An integer to specify the number of posterior sampling iteration of each chain in MCMC for stan models. Default is 1000. This parameter only works for Bayesian models.
#' @param nchain An integer to specify the number of MCMC chains for stan models. Default is 4. This parameter only works for Bayesian models.
#' @param HPD_interval A logical scalar; if true, the calculated credible intervals for stan models are highest posterior density intervals. Otherwise the intervals are symmetric. Default is false. This parameter only works for Bayesian models.
#'
#' @rawNamespace import(stats, except = filter)
#' @import rstanarm
#' @import survey
#' @import mgcv
#' @importFrom dplyr mutate_at filter %>% pull
#' @import stringr
#' @import BART
#'
#'
#' @return A list. Each element contains the estimate and CIs for a subset or the whole data analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' ## simulate data from the 'simulate' function, with nonlinear association setting 3.
#' ## The continuous variable X is discretized into categorical variable auX_3 with 3 categories,
#' ## and auX_10 with 10 categories.
#' data = simulate(N = 3000, discretize = c(3, 10), setting = 3, seed = 123)
#'
#' ## The dataset consists of a continuous outcome Y, 3 binary variables Z1, Z2, Z3;
#' ## discretized variables auX_3 and auX_10.
#' ## Propensity scores true_pi and it log transformation logit_true_pi are calculated.
#' population = data$population # get population, 3000 cases
#' samples = data$samples # get samples, about 600 cases
#' ipw = 1 / samples$true_pi # get the inverse probability weighting
#'
#' ## True value
#' true_mean = mean(population$Y1)
#'
#' ## IPW sample mean, with analysis on subset Z1 == 1 & Z2 == 1.
#' IPW_sample_mean = auxsurvey("~Y1",  auxiliary = NULL, weights = ipw,
#'                    samples = samples, population = population,
#'                    subset = c("Z1 == 1 & Z2 == 1"), method = "sample_mean", levels = 0.95)
#'
#' ## rake, with analysis on subsets Z1 == 1 and Z1 == 1 & Z2 == 1.
#' rake = auxsurvey("~Y1",  auxiliary = "Z1 + Z2 + Z3 + auX_10", samples = samples,
#'          population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1"),
#'          method = "rake", levels = 0.95)
#'
#' ## IPW post-stratification, no subset analysis.
#' IPW_postStratify3 = auxsurvey("~Y1",  auxiliary = "Z1 + Z2 + Z3 + auX_3", weights = ipw,
#'                      samples = samples, population = population,
#'                      method = "postStratify", levels = 0.95)
#'
#' ## MRP, with analysis on subsets Z1 == 1, Z1 == 1 & Z2 == 1, Z1 == 1 & Z2 == 1 & Z3 == 1.
#' MRP = auxsurvey("Y1~1 + Z1",  auxiliary = "Z2 + Z3:auX_10", samples = samples,
#'         population = population,
#'         subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"),
#'         method = "MRP", levels = 0.95, nskip = 4000, npost = 4000,
#'         nchain = 1, stan_verbose = F, HPD_interval = T)
#'
#' ## GAMP, no subset analysis.
#' GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxiliary = "s(auX_10) + s(logit_true_pi, by = Z1)",
#'           samples = samples, population = population, subset = NULL, method = "GAMP",
#'          levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
#'
#'
#' ## BART with analysis on subsets Z1 == 1, Z1 == 1 & Z2 == 1, Z1 == 1 & Z2 == 1 & Z3 == 1.
#' MRP = auxsurvey("Y1~1 + Z1 + Z2 + Z3 + auX_10",  auxiliary = NULL, samples = samples,
#'         population = population,
#'         subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"),
#'         method = "BART", levels = 0.95, nskip = 4000, npost = 4000,
#'         nchain = 1, stan_verbose = F, HPD_interval = T)
#' }
auxsurvey <- function(formula, auxiliary = NULL, samples, population = NULL, subset = NULL, family = gaussian(), method = c("sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear", "BART"), weights = NULL, levels = c(0.95, 0.8, 0.5), stan_verbose = TRUE, show_plot = TRUE, nskip = 1000, npost = 1000, nchain = 4, HPD_interval = FALSE){
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

  #auxiliary = str_trim(str_split(auxiliary, "\\+", simplify = T))
  #auxiliary = str_remove_all()
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
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    auxiliary = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T))
    auxiliary = auxiliary[!is.na(auxiliary) & auxiliary != ""]
    return(rake_wt(samples, population, auxiliary, svyVar, subset, family = family, levels, weights, maxiter = 50))
  }
  if(method == "postStratify"){
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
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
    covariates = stringr::str_trim(stringr::str_split(stringr::str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    samples = dplyr::mutate_at(samples, all.vars(as.formula(auxiliary)), as.factor)
    population = dplyr::mutate_at(population, all.vars(as.formula(auxiliary)), as.factor)
    auxiliary = stringr::str_replace_all(auxiliary, "\\*", ":")
    auxiliary = stringr::str_split_i(as.character(auxiliary), "~", 2)
    auxiliary = stringr::str_split(auxiliary, "\\+", simplify = T)

    if(length(covariates) == 0){
      outcome_formula = paste0(paste0(svyVar, "~"), paste0("(1|", auxiliary, ")", collapse = "+"))
    }else{
      outcome_formula = paste(paste0(svyVar, "~", paste0(covariates, collapse = "+")), paste0("(1|", auxiliary, ")", collapse = "+"), sep= "+")
    }
    cat("The formula for the MRP model is ", outcome_formula, "\n")
    MRP_est = svyBayesmod(samples, population, outcome_formula, "stan_glmer", subset, family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = show_plot, useTrueSample = F, stan_verbose = stan_verbose, shortest_CI = HPD_interval)
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
    #covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    #if("." %in% covariates){
    #  covariates = setdiff(names(samples), svyVar)
    #}

    #samples_matrix = model.matrix(as.formula(paste(paste0("~", str_split_i(formula, "~", 2)), paste0(auxiliary, collapse = "+"), sep = "+")), data = samples)
    samples = mutate_at(samples, all.vars(as.formula(auxiliary)), as.numeric)
    population = mutate_at(population, intersect(all.vars(as.formula(auxiliary)), colnames(population)), as.numeric)

    auxiliary = stringr::str_replace_all(auxiliary, "\\*", ":")
    auxiliary_fixed = setdiff(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxiliary)))[stringr::str_detect(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), "s\\(.*\\)")]


    outcome_formula = paste(formula, paste(auxiliary_fixed, collapse = "+"), sep = "+")
    auxiliary_random =  union(intersect(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxiliary))), stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T)[is.na(stringr::str_match(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), "s\\(.*\\)"))])

    if(length(auxiliary_random) > 0){
      samples = mutate_at(samples, all.vars(as.formula(paste0("~", auxiliary_random))), as.factor)
      population = mutate_at(population, all.vars(as.formula(paste0("~", auxiliary_random))), as.factor)
    }
    if(length(auxiliary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxiliary_random, ")", collapse = "+")))
      cat("The formula for the GAMP model is ", paste0(outcome_formula[1], str_replace(outcome_formula[2], "~", "+"), collapse = ""), "\n")
    }else{
      outcome_formula = c(outcome_formula, NULL)
      cat("The formula for the GAMP model is ", outcome_formula[1], "\n")
    }

    #outcome_formula = paste(formula, paste0(str_remove_all(auxiliary, "~"), collapse = "+"), paste0("(1|", auxiliary_random, ")", collapse = "+"),  sep = "+")

    GAMP_est = svyBayesmod(samples, population, outcome_formula, "stan_gamm4", subset, family, levels, weights, nskip, npost, nchain, printmod = T, doFigure = F, useTrueSample = T, stan_verbose = stan_verbose, shortest_CI = HPD_interval)
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
    population = mutate_at(population, intersect(all.vars(as.formula(auxiliary)), colnames(population)), as.numeric)


    auxiliary_fixed = setdiff(stringr::str_split(stringr::str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxiliary)))


    outcome_formula = paste(formula, paste(auxiliary_fixed, collapse = "+"), collapse = "+")
    auxiliary_random = intersect(str_split(str_split_i(as.character(auxiliary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxiliary)))
    samples = mutate_at(samples, auxiliary_random, as.factor)
    population = mutate_at(population, auxiliary_random, as.factor)
    if(length(auxiliary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxiliary_random, ")", collapse = "+")))
      cat("The formula for the linear model is ", paste0(outcome_formula[1], str_replace(outcome_formula[2], "~", "+"), collapse = ""), "\n")
      Linear_est = svyBayesmod(samples, population, outcome_formula, "stan_glmer", family, levels, weights, nskip, npost, nchain, printmod = T, doFigure = F, useTrueSample = T, stan_verbose = stan_verbose, shortest_CI = HPD_interval)
    }else{
      outcome_formula = c(outcome_formula, NULL)
      cat("The formula for the linear model is ", outcome_formula, "\n")
      Linear_est = svyBayesmod(samples, population, outcome_formula, "stan_glm", subset, family, levels, weights, nskip, npost, nchain, printmod = T, doFigure = F, useTrueSample = T, stan_verbose = stan_verbose, shortest_CI = HPD_interval)
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
    covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
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
      X_train = as.matrix(samples[, covariates])
      y_train = dplyr::pull(samples, svyVar)
      if(!stan_verbose)
        model <- BART::pbart(X_train, y_train, ndpost = npost, nskip = nskip)
      else{
        model <- BART::pbart(X_train, y_train, ndpost = npost, nskip = nskip)
      }
    }
    if(family$family == "gaussian"){
      X_train = as.matrix(samples[, covariates])
      y_train = dplyr::pull(samples, svyVar)
      if(!stan_verbose)
        model <- BART::wbart(X_train, y_train, ndpost = npost, nskip = nskip)
      else{
        model <- BART::wbart(X_train, y_train, ndpost = npost, nskip = nskip)
      }
    }
    subset = c("T", subset)
    infr = sapply(subset, function(s){
      svypopu1 = as.matrix(dplyr::filter(population, eval(parse(text = s)))[, covariates])
      svysmpl1 = as.matrix(dplyr::filter(samples, eval(parse(text = s)))[, covariates])

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
          ci = coda::HPDinterval(post_est, prob = level, names = T)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1))
          ci = quantile(post_est, probs = ci, names = T)
          names(ci) = str_replace(names(ci), "%", " %")
        }
        ci
      }, simplify = F)

      tCI = do.call("c", tCI)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(svysmpl1), population_size = nrow(svypopu1)))
      return(infr)
    }, simplify = F)
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

