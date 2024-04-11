library(gtools)
library(scales)
library(coda)
library(survey)
library(pROC)
library(tidyverse)
library(mgcv)
library(rstanarm)
library(rlang)

uwt <- function(svysmpl, svyVar, subset = NULL, invlvls, weights = NULL) {
  subset = c("T", subset)
  infr = sapply(subset, function(s){


    index <- eval(rlang::parse_expr(s), envir = svysmpl)
    a = svysmpl[index, ]
    weights1 = weights[index]
    obj <- paste(svyVar, 1, sep = '~') %>% as.formula(.) %>%
      lm(., data = a, weights = weights1)
    invls <- sapply(invlvls, function(lvl){
      inv <- confint(object = obj, level = lvl);
      rownames(inv) <- lvl;
      return(inv)
    }, simplify = F)
    invls <- invls %>% do.call('cbind', .)
    infr <- cbind(est = obj$coefficients, se = sqrt(diag(vcov(obj))), invls, sample_size = length(obj$fitted.values))
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

rake_wt <- function(svysmpl, svypopu, auxVars, svyVar, subset = NULL, family = gaussian(), invlvls, weights = NULL, maxiter = 50) {
  subset = c("T", subset)
  svysmpl$fpc <- nrow(svypopu)
  if(is.null(weights))
    des <- svydesign(id = ~1, weights = ~1, data = svysmpl, fpc = ~fpc)
  else{
    des <- svydesign(id = ~1, weights = ~weights, data = svysmpl, fpc = ~fpc)
  }
  popu_tab <- lapply(auxVars, function (Var) {
    Var = paste(all.vars(as.formula(paste0("~", Var))), collapse = "+")
    fmla <- paste('~', Var, sep = ' ') %>% as.formula(.)
    tab <- xtabs(fmla, svypopu)
    return(tab)
  } )
  fmla_list <- lapply(auxVars, function (Var) {
    Var = paste(all.vars(as.formula(paste0("~", Var))), collapse = "+")
    paste('~', Var, sep = ' ') %>% as.formula(.) %>% return(.)
  } )
  rakingobj <- rake(des, fmla_list, popu_tab,
                    control = list(maxit = maxiter, epsilon = 1, verbose = FALSE) )
  infr = sapply(subset, function(s){
    rakingobj = subset(rakingobj, eval(parse(text = s)))
    if(family$family == "binomial"){
      suppressWarnings(rakest <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula(.) %>%
        svyciprop(., rakingobj, method = "logit", level = lv), simplify = F))
      tCI <- lapply(rakest, function(i){
        ci = confint(i, df = degf(rakingobj), parm = svyVar)
        colnames(ci) = str_replace(colnames(ci), "%", " %")
        ci
      })%>% do.call("cbind", .)
      rakest = rakest[[1]]
    }
    if(family$family == "gaussian"){
      rakest <- paste0('~', svyVar) %>% as.formula(.) %>%
        svymean(., rakingobj)
      tCI <- sapply(invlvls, confint, object = rakest, df = degf(rakingobj),
                    parm = svyVar, simplify = F) %>% do.call("cbind", .)
    }
    infr <- cbind(est = rakest[svyVar], se = sqrt(diag(vcov(rakest))), tCI, sample_size = degf(rakingobj) + 1, population_size = nrow(filter(svypopu, eval(parse(text = s)))))
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

postStr_wt <- function(svysmpl, svypopu, auxVars, svyVar, subset = NULL, family = gaussian(), invlvls, weights = NULL) {
# svysmpl is sample
# svypopu is population
# auxVar is a vector of covariates names for post-stratification
# svyVar is the outcome variables
# invlvls is a vector of confidence levels
  subset = c("T", subset)
  svysmpl$fpc <- nrow(svypopu) # get population sample size

  # specify svydesign
  if(is.null(weights))
    des <- svydesign(id = ~1, weights = ~1, data = svysmpl, fpc = ~fpc)
  else{
    des <- svydesign(id = ~1, weights = ~weights, data = svysmpl, fpc = ~fpc)
  }


  fmla <- paste(auxVars, collapse = "+") %>% paste('~', ., sep = ' ') %>% as.formula(.)

  # set for post-stratification
  tab <- xtabs(fmla, svypopu)
  PSobj <- postStratify(des, fmla, tab, partial = TRUE)
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
      suppressWarnings(PSest <- sapply(invlvls, function(lv) paste0('~', svyVar) %>% as.formula(.) %>%
                         svyciprop(., PSobj, method = "logit", level = lv), simplify = F))
      tCI <- lapply(PSest, function(i){
        ci = confint(i, df = degf(PSobj), parm = svyVar)
        colnames(ci) = str_replace(colnames(ci), "%", " %")
        ci
      })%>% do.call("cbind", .)
      PSest = PSest[[1]]
    }
    if(family$family == "gaussian"){
      PSest <- paste0('~', svyVar) %>% as.formula(.) %>% svymean(., PSobj)
      tCI <- sapply(invlvls, confint, object = PSest, df = degf(PSobj),
                    parm = svyVar, simplify = F) %>% do.call("cbind", .)
    }

    # get estmates and standard error
    infr <- cbind(est = PSest[svyVar], se = sqrt(diag(vcov(PSest))), tCI, sample_size = degf(PSobj) + 1, population_size = nrow(filter(svypopu, eval(parse(text = s)))))
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


svyBayesmod <- function(svysmpl, svypopu, outcome_formula, BayesFun, subset = NULL, family = gaussian(), invlvls, weights = NULL, nskip = 1000, npost = 1000, nchain = 4, printmod = TRUE, doFigure = FALSE, useTrueSample = F, stan_verbose = F, shortest_CI = F) {
  #print(outcome_formula)
  subset = c("T", subset)
  fmla <- outcome_formula[1] %>% as.formula(.)
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

  if (printmod) summary(bayesmod) %>% print(.)

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
      svypopu1 = filter(svypopu, eval(parse(text = s)))
      yhats <- posterior_epred(bayesmod, svypopu1)

      post_est = rowMeans(yhats)

      tCI = sapply(invlvls, function(level){
        if(shortest_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, probs = level, names = T)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1)) %>%
            quantile(post_est, probs = ., names = T)
          names(ci) = str_replace(names(ci), "%", " %")
        }

        ci
      }, simplify = F) %>% do.call("c", .)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(filter(svysmpl, eval(parse(text = s)))), population_size = nrow(svypopu1)))
    }, simplify = F)
    names(infr)[1] = "All"
    if(length(infr) == 1)
      return(infr[[1]])
    return(infr)
  }else{
    infr = sapply(subset, function(s){
      svypopu1 = filter(svypopu, eval(parse(text = s)))
      svysmpl1 = filter(svysmpl, eval(parse(text = s)))

      yhats <- posterior_predict(bayesmod, setdiff(svypopu1, svysmpl1), re.form = NA) # npost, non-sample_size

      yhat_tot <- yhats %>% apply(., 1, sum) # npost * 1
      yobs_tot <- sum(svysmpl1 %>% dplyr::select(str_split_1(outcome_formula[1], "~")[1])) # a single value
      post_est = (yhat_tot + yobs_tot) / nrow(svypopu1) #npost
      tCI = sapply(invlvls, function(level){
        #confint(post_est, level = 0.95)
        if(shortest_CI){
          class(post_est) <- 'mcmc'
          ci = coda::HPDinterval(post_est, prob = level, names = T)
          names(ci) = paste(((1 - level)/2 * c(1, -1) + c(0, 1)) * 100, "%")
        }else{
          ci = ((1 - level)/2 * c(1, -1) + c(0, 1)) %>%
            quantile(post_est, probs = ., names = T)
          names(ci) = str_replace(names(ci), "%", " %")
        }

        ci
      }, simplify = F) %>% do.call("c", .)
      infr <- rbind(c(post_mean_est = mean(post_est), post_median_est = median(post_est), se = sd(post_est), tCI, sample_size = nrow(svysmpl1), population_size = nrow(svypopu1)))
      return(infr)
    }, simplify = F)
    names(infr)[1] = "All"
    return(infr)
  }
}



auxsurvey <- function(formula, auxilary = NULL, samples, population = NULL, subset = NULL, family = gaussian(), method = c("sample_mean", "rake", "postStratify", "MRP", "GAMP", "linear"), weights = NULL, levels = c(0.95, 0.8, 0.5), stan_verbose = T, show_plot = T, nskip = 1000, npost = 1000, nchain = 4, HPD_interval = F, ...){
  svyVar = str_trim(str_split_1(as.character(formula), "~"))
  svyVar = svyVar[svyVar != ""][1]
  if(!is.null(auxilary)){
    auxilary = as.character(auxilary)
    auxilary = str_remove_all(auxilary, " ")
    auxilary = paste0(auxilary, collapse = "+")
    if(!str_detect(auxilary, "^~")){
      auxilary = paste0("~",auxilary)
    }
  }
  if(length(setdiff(union(all.vars(as.formula(auxilary)), all.vars(as.formula(formula))), colnames(population))) > 0){
    stop(paste0("unidentified variables: ", paste0(setdiff(union(all.vars(as.formula(auxilary)), all.vars(as.formula(formula))), colnames(population)), collapse = ", "), collapse = ", "))
  }

  #auxilary = str_trim(str_split(auxilary, "\\+", simplify = T))
  #auxilary = str_remove_all()
  method = match.arg(method)
  # if(is.null(propensity_score))
  #   weights = NULL
  # else{
  #   weights = 1 / propensity_score
  # }
  if(method == "sample_mean"){
    return(uwt(samples, svyVar, subset, levels, weights))
  }
  if(method == "rake"){
    covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    auxilary = str_trim(str_split(str_split_i(as.character(auxilary), "~", 2), "\\+", simplify = T))
    auxilary = auxilary[!is.na(auxilary) & auxilary != ""]
    return(rake_wt(samples, population, auxilary, svyVar, subset, family = family, levels, weights, maxiter = 50))
  }
  if(method == "postStratify"){
    covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
      covariates = setdiff(names(samples), svyVar)
    }
    auxilary = unique(c(all.vars(as.formula(auxilary)), covariates))
    auxilary = auxilary[!is.na(auxilary) & auxilary != ""]
    return(postStr_wt(samples, population, auxilary, svyVar, subset, family = family, levels, weights))
  }
  if(method == "MRP"){
    if(is.null(nskip)) nskip = 1000
    if(is.null(npost)) npost = 1000
    if(is.null(nchain)) nchain = 4
    covariates = str_trim(str_split(str_split_i(as.character(formula), "~", 2), "\\+", simplify = T))
    covariates = setdiff(covariates, svyVar)
    if("." %in% covariates){
     covariates = setdiff(names(samples), svyVar)
    }
    samples = mutate_at(samples, all.vars(as.formula(auxilary)), as.factor)
    population = mutate_at(population, all.vars(as.formula(auxilary)), as.factor)
    auxilary = str_replace_all(auxilary, "\\*", ":")
    auxilary = str_split_i(as.character(auxilary), "~", 2)
    auxilary = str_split(auxilary, "\\+", simplify = T)

    if(length(covariates) == 0){
      outcome_formula = paste0(paste0(svyVar, "~"), paste0("(1|", auxilary, ")", collapse = "+"))
    }else{
      outcome_formula = paste(paste0(svyVar, "~", paste0(covariates, collapse = "+")), paste0("(1|", auxilary, ")", collapse = "+"), sep= "+")
    }
    cat("The formula for the MRP model is ", outcome_formula, "\n")
    MRP_est = svyBayesmod(samples, population, outcome_formula, "stan_glmer", subset, family, levels, weights, nskip, npost, nchain, printmod = TRUE, doFigure = show_plot, useTrueSample = F, stan_verbose = stan_verbose, shortest_CI = HPD_interval)
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

    #samples_matrix = model.matrix(as.formula(paste(paste0("~", str_split_i(formula, "~", 2)), paste0(auxilary, collapse = "+"), sep = "+")), data = samples)
    samples = mutate_at(samples, all.vars(as.formula(auxilary)), as.numeric)
    population = mutate_at(population, intersect(all.vars(as.formula(auxilary)), colnames(population)), as.numeric)


    auxilary_fixed = setdiff(str_split(str_split_i(as.character(auxilary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxilary)))


    outcome_formula = paste(formula, paste(auxilary_fixed, collapse = "+"), sep = "+")
    auxilary_random = intersect(str_split(str_split_i(as.character(auxilary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxilary)))
    samples = mutate_at(samples, auxilary_random, as.factor)
    population = mutate_at(population, auxilary_random, as.factor)
    if(length(auxilary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxilary_random, ")", collapse = "+")))
      cat("The formula for the GAMP model is ", paste0(outcome_formula[1], str_replace(outcome_formula[2], "~", "+"), collapse = ""), "\n")
    }else{
      outcome_formula = c(outcome_formula, NULL)
      cat("The formula for the GAMP model is ", outcome_formula[1], "\n")
    }

    #outcome_formula = paste(formula, paste0(str_remove_all(auxilary, "~"), collapse = "+"), paste0("(1|", auxilary_random, ")", collapse = "+"),  sep = "+")

    GAMP_est = svyBayesmod(samples, population, outcome_formula, "stan_gamm4", subset, family, levels, weights, nskip, npost, nchain, printmod = T, doFigure = F, useTrueSample = T, stan_verbose = stan_verbose, shortest_CI = HPD_interval)


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
    samples = mutate_at(samples, all.vars(as.formula(auxilary)), as.numeric)
    population = mutate_at(population, intersect(all.vars(as.formula(auxilary)), colnames(population)), as.numeric)


    auxilary_fixed = setdiff(str_split(str_split_i(as.character(auxilary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxilary)))


    outcome_formula = paste(formula, paste(auxilary_fixed, collapse = "+"), collapse = "+")
    auxilary_random = intersect(str_split(str_split_i(as.character(auxilary), "~", 2), "\\+", simplify = T), all.vars(as.formula(auxilary)))
    samples = mutate_at(samples, auxilary_random, as.factor)
    population = mutate_at(population, auxilary_random, as.factor)
    if(length(auxilary_random) != 0){
      outcome_formula = c(outcome_formula, paste0("~", paste0("(1|", auxilary_random, ")", collapse = "+")))
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

    #outcome_formula = paste(formula, paste0(str_remove_all(auxilary, "~"), collapse = "+"), paste0("(1|", auxilary_random, ")", collapse = "+"),  sep = "+")
    return(Linear_est)
  }
}

