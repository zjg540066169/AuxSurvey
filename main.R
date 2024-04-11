library(plyr)
library(tidyverse)
source("./R/simulate.R")
source("./R/AuxSurvey.R")

args = commandArgs(trailingOnly=TRUE)
seed = as.integer(args[1])
setting = as.integer(args[2])

seed = 1
setting = 1

library(doParallel)
library(foreach)
registerDoParallel(6)

setting = 2
a2 = foreach(seed = 1:500, .combine = c, .verbose = T) %dopar%{
  set.seed(seed)
  data = simulate(N = 3000, n = 600, discretize = c(3, 5, 10), setting = setting, seed = seed)
  population = data$population
  samples = data$samples



  true_theta_1 = mean(population$Y1)
  true_theta_1_Z1_1 = mean(unlist(population[population$Z1 == 1, "Y1"]))
  true_theta_1_Z1_1_Z2_1 = mean(unlist(population[population$Z1 == 1 & population$Z2 == 1, "Y1"]))
  true_theta_1_Z1_1_Z2_1_Z3_1 = mean(unlist(population[population$Z1 == 1 & population$Z2 == 1 & population$Z3 == 1, "Y1"]))

  true_theta_2 = mean(population$Y2)
  true_theta_2_Z1_1 = mean(unlist(population[population$Z1 == 1, "Y2"]))
  true_theta_2_Z1_1_Z2_1 = mean(unlist(population[population$Z1 == 1 & population$Z2 == 1, "Y2"]))
  true_theta_2_Z1_1_Z2_1_Z3_1 = mean(unlist(population[population$Z1 == 1 & population$Z2 == 1 & population$Z3 == 1, "Y2"]))

  if(setting == 1){
    # results for Y1
    continuous = list(
      sample_mean = auxsurvey("~Y1",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      MRP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(logit_true_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(logit_estimated_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y1~1 + Z1 + Z2 + Z3 + X + I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )
    # results for Y2
    binary = list(
      sample_mean = auxsurvey("~Y2",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      MRP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(logit_true_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(logit_estimated_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y2~1 + Z1 + Z2 + Z3 + X + I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )


  }

  if(setting == 2){
    # results for Y1
    continuous = list(
      sample_mean = auxsurvey("~Y1",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake10 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      MRP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "auX_10 + auW_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(auW_10) + s(logit_true_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(auW_10) + s(logit_estimated_pi_with_W)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y1~1 + Z1 + Z2 + Z3 + X + I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )

    # results for Y2
    binary = list(
      sample_mean = auxsurvey("~Y2",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5 + auW_5", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake10 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_10 + auW_10", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3 + auW_3", weights = 1 / samples$estimated_pi_with_W, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      MRP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "auX_10 + auW_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(auW_10) + s(logit_true_pi)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10) + s(auW_10) + s(logit_estimated_pi_with_W)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y2~1 + Z1 + Z2 + Z3 + X + I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )

  }

  if(setting == 3){
    # results for Y1
    continuous = list(
      sample_mean = auxsurvey("~Y1",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y1",  auxilary = NULL, weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      rake5_2 = auxsurvey("~Y1",  auxilary = "Z2 + Z3 + Z1 * auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      true_IPW_rake5_2 = auxsurvey("~Y1",  auxilary = "Z2 + Z3 + Z1 * auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y1", auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),
      estimated_IPW_rake5_2 = auxsurvey("~Y1",  auxilary = "Z2 + Z3 + Z1 * auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "rake", levels = 0.95),

      postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y1",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "postStratify", levels = 0.95),
      MRP_interaction = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "Z1 * auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      MRP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10, by = Z1) + s(logit_true_pi, by = Z1)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y1~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10, by = Z1) + s(logit_estimated_pi, by = Z1)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y1~1 + Z1 + Z2 + Z3 + X + I(X^2) + Z1 * X + Z1 * I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )

    # results for Y2
    binary = list(
      sample_mean = auxsurvey("~Y2",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      true_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      estimated_IPW_sample_mean = auxsurvey("~Y2",  auxilary = NULL, weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "sample_mean", levels = 0.95),
      rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake5 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      rake5_2 = auxsurvey("~Y2",  auxilary = "Z2 + Z3 + Z1 * auX_5", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      true_IPW_rake5_2 = auxsurvey("~Y2",  auxilary = "Z2 + Z3 + Z1 * auX_5", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake5 = auxsurvey("~Y2", auxilary = "Z1 + Z2 + Z3 + auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      estimated_IPW_rake5_2 = auxsurvey("~Y2",  auxilary = "Z2 + Z3 + Z1 * auX_5", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "rake", levels = 0.95),
      postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      true_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$true_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      estimated_IPW_postStratify3 = auxsurvey("~Y2",  auxilary = "Z1 + Z2 + Z3 + auX_3", weights = 1 / samples$estimated_pi, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "postStratify", levels = 0.95),
      MRP_interaction = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "Z1 * auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      MRP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "auX_10", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "MRP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10, by = Z1) + s(logit_true_pi, by = Z1)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      estimated_GAMP = auxsurvey("Y2~1 + Z1 + Z2 + Z3",  auxilary = "s(auX_10, by = Z1) + s(logit_estimated_pi, by = Z1)", samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "GAMP", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T),
      true_model = auxsurvey("Y2~1 + Z1 + Z2 + Z3 + X + I(X^2) + Z1 * X + Z1 * I(X^2)",  auxilary = NULL, samples = samples, population = population, subset = c("Z1 == 1", "Z1 == 1 & Z2 == 1", "Z1 == 1 & Z2 == 1 & Z3 == 1"), family = binomial(link = "logit"), method = "linear", levels = 0.95, nskip = 4000, npost = 4000, nchain = 1, stan_verbose = F, HPD_interval = T)
    )

  }

  continuous = lapply(continuous, function(result){
    result[["All"]] = result[["All"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_1 <= `97.5 %` & true_theta_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_1, est - true_theta_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_1
      )

    result[["Z1 == 1"]] = result[["Z1 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_1_Z1_1 <= `97.5 %` & true_theta_1_Z1_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_1_Z1_1, est - true_theta_1_Z1_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_1_Z1_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_1_Z1_1
      )

    result[["Z1 == 1 & Z2 == 1"]] = result[["Z1 == 1 & Z2 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_1_Z1_1_Z2_1 <= `97.5 %` & true_theta_1_Z1_1_Z2_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_1_Z1_1_Z2_1, est - true_theta_1_Z1_1_Z2_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_1_Z1_1_Z2_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_1_Z1_1_Z2_1
      )

    result[["Z1 == 1 & Z2 == 1 & Z3 == 1"]] = result[["Z1 == 1 & Z2 == 1 & Z3 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_1_Z1_1_Z2_1_Z3_1 <= `97.5 %` & true_theta_1_Z1_1_Z2_1_Z3_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_1_Z1_1_Z2_1_Z3_1, est - true_theta_1_Z1_1_Z2_1_Z3_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_1_Z1_1_Z2_1_Z3_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_1_Z1_1_Z2_1_Z3_1
      )
    result
  })
  binary = lapply(binary, function(result){
    #print(names(result))
    result[["All"]] = result[["All"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_2 <= `97.5 %` & true_theta_2 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_2, est - true_theta_2
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_2, NA
             ),
             square_error = bias^2,
             theta = true_theta_2
      )

    result[["Z1 == 1"]] = result[["Z1 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_2_Z1_1 <= `97.5 %` & true_theta_2_Z1_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_2_Z1_1, est - true_theta_2_Z1_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_2_Z1_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_2_Z1_1
      )

    result[["Z1 == 1 & Z2 == 1"]] = result[["Z1 == 1 & Z2 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_2_Z1_1_Z2_1 <= `97.5 %` & true_theta_2_Z1_1_Z2_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_2_Z1_1_Z2_1, est - true_theta_2_Z1_1_Z2_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_2_Z1_1_Z2_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_2_Z1_1_Z2_1
      )
    result[["Z1 == 1 & Z2 == 1 & Z3 == 1"]] = result[["Z1 == 1 & Z2 == 1 & Z3 == 1"]] %>%
      as_tibble() %>%
      mutate(AW = `97.5 %` - `2.5 %`,
             coverage = ifelse(
               true_theta_2_Z1_1_Z2_1_Z3_1 <= `97.5 %` & true_theta_2_Z1_1_Z2_1_Z3_1 >= `2.5 %`, T, F
             ),
             bias = ifelse(
               !"est" %in% colnames(.), post_median_est - true_theta_2_Z1_1_Z2_1_Z3_1, est - true_theta_2_Z1_1_Z2_1_Z3_1
             ),
             bias_mean = ifelse(
               !"est" %in% colnames(.), post_mean_est - true_theta_2_Z1_1_Z2_1_Z3_1, NA
             ),
             square_error = bias^2,
             theta = true_theta_2_Z1_1_Z2_1_Z3_1
      )
    result
  })


  model_names = names(continuous)
  continuous_All = do.call(rbind.fill, lapply(continuous, function(c) c$All))
  rownames(continuous_All) = model_names
  continuous_Z1_1 = do.call(rbind.fill, lapply(continuous, function(c) c[["Z1 == 1"]]))
  rownames(continuous_Z1_1) = model_names
  continuous_Z1_1_Z2_1 = do.call(rbind.fill, lapply(continuous, function(c) c[["Z1 == 1 & Z2 == 1"]]))
  rownames(continuous_Z1_1_Z2_1) = model_names
  continuous_Z1_1_Z2_1_Z3_1 = do.call(rbind.fill, lapply(continuous, function(c) c[["Z1 == 1 & Z2 == 1 & Z3 == 1"]]))
  rownames(continuous_Z1_1_Z2_1_Z3_1) = model_names
  model_names = names(binary)
  binary_All = do.call(rbind.fill, lapply(binary, function(c) c$All))
  rownames(binary_All) = model_names
  binary_Z1_1 = do.call(rbind.fill, lapply(binary, function(c) c[["Z1 == 1"]]))
  rownames(binary_Z1_1) = model_names
  binary_Z1_1_Z2_1 = do.call(rbind.fill, lapply(binary, function(c) c[["Z1 == 1 & Z2 == 1"]]))
  rownames(binary_Z1_1_Z2_1) = model_names
  binary_Z1_1_Z2_1_Z3_1 = do.call(rbind.fill, lapply(binary, function(c) c[["Z1 == 1 & Z2 == 1 & Z3 == 1"]]))
  rownames(binary_Z1_1_Z2_1_Z3_1) = model_names
  result = list(list(continuous_All = continuous_All,
                              continuous_Z1_1 = continuous_Z1_1,
                              continuous_Z1_1_Z2_1 = continuous_Z1_1_Z2_1,
                              continuous_Z1_1_Z2_1_Z3_1 = continuous_Z1_1_Z2_1_Z3_1,
                              binary_All = binary_All,
                              binary_Z1_1 = binary_Z1_1,
                              binary_Z1_1_Z2_1 = binary_Z1_1_Z2_1,
                              binary_Z1_1_Z2_1_Z3_1 = binary_Z1_1_Z2_1_Z3_1
  ))
  names(result) = seed
  return(result)
}



summary_results = function(a3, setting = 1){
  result_ALL = do.call("rbind", sapply(1:length(a3), function(i){
    a3[[i]]$continuous_All$replicates = i
    a3[[i]]$continuous_All$outcome = "continuous"
    a3[[i]]$continuous_All = rownames_to_column(a3[[i]]$continuous_All, 'model')

    a3[[i]]$binary_All$replicates = i
    a3[[i]]$binary_All$outcome = "binary"
    a3[[i]]$binary_All = rownames_to_column(a3[[i]]$binary_All, 'model')

    rbind(a3[[i]]$continuous_All, a3[[i]]$binary_All)
  }, simplify = F)) %>%
    mutate(
      model = factor(model, levels = rownames(a3[[1]]$continuous_All)),
      outcome = factor(outcome, levels = c("continuous", "binary"))
    )


  result_Z1_1 = do.call("rbind", sapply(1:length(a3), function(i){
    a3[[i]]$continuous_Z1_1$replicates = i
    a3[[i]]$continuous_Z1_1$outcome = "continuous"
    a3[[i]]$continuous_Z1_1 = rownames_to_column(a3[[i]]$continuous_Z1_1, 'model')

    a3[[i]]$binary_Z1_1$replicates = i
    a3[[i]]$binary_Z1_1$outcome = "binary"
    a3[[i]]$binary_Z1_1 = rownames_to_column(a3[[i]]$binary_Z1_1, 'model')

    rbind(a3[[i]]$continuous_Z1_1, a3[[i]]$binary_Z1_1)
  }, simplify = F)) %>%
    mutate(
      model = factor(model, levels = rownames(a3[[1]]$continuous_Z1_1)),
      outcome = factor(outcome, levels = c("continuous", "binary"))
    )

  result_Z1_1_Z2_1 = do.call("rbind", sapply(1:length(a3), function(i){
    a3[[i]]$continuous_Z1_1_Z2_1$replicates = i
    a3[[i]]$continuous_Z1_1_Z2_1$outcome = "continuous"
    a3[[i]]$continuous_Z1_1_Z2_1 = rownames_to_column(a3[[i]]$continuous_Z1_1_Z2_1, 'model')

    a3[[i]]$binary_Z1_1_Z2_1$replicates = i
    a3[[i]]$binary_Z1_1_Z2_1$outcome = "binary"
    a3[[i]]$binary_Z1_1_Z2_1 = rownames_to_column(a3[[i]]$binary_Z1_1_Z2_1, 'model')

    rbind(a3[[i]]$continuous_Z1_1_Z2_1, a3[[i]]$binary_Z1_1_Z2_1)
  }, simplify = F)) %>%
    mutate(
      model = factor(model, levels = rownames(a3[[1]]$continuous_Z1_1_Z2_1)),
      outcome = factor(outcome, levels = c("continuous", "binary"))
    )

  result_Z1_1_Z2_1_Z3_1 = do.call("rbind", sapply(1:length(a3), function(i){
    a3[[i]]$continuous_Z1_1_Z2_1_Z3_1$replicates = i
    a3[[i]]$continuous_Z1_1_Z2_1_Z3_1$outcome = "continuous"
    a3[[i]]$continuous_Z1_1_Z2_1_Z3_1 = rownames_to_column(a3[[i]]$continuous_Z1_1_Z2_1_Z3_1, 'model')

    a3[[i]]$binary_Z1_1_Z2_1_Z3_1$replicates = i
    a3[[i]]$binary_Z1_1_Z2_1_Z3_1$outcome = "binary"
    a3[[i]]$binary_Z1_1_Z2_1_Z3_1 = rownames_to_column(a3[[i]]$binary_Z1_1_Z2_1_Z3_1, 'model')

    rbind(a3[[i]]$continuous_Z1_1_Z2_1_Z3_1, a3[[i]]$binary_Z1_1_Z2_1_Z3_1)
  }, simplify = F)) %>%
    mutate(
      model = factor(model, levels = rownames(a3[[1]]$continuous_Z1_1_Z2_1_Z3_1)),
      outcome = factor(outcome, levels = c("continuous", "binary"))
    )


  print(result_ALL %>%
    ggplot() + geom_boxplot(aes(x = model, y = bias)) + geom_hline(yintercept=0, linetype="dashed", color = "red") + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Bias"))
  )


  print(result_ALL %>%
    ggplot() + geom_boxplot(aes(x = model, y = square_error)) + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Square Error"))
  )

  print(result_ALL %>%
    ggplot() + geom_boxplot(aes(x = model, y = AW)) + facet_grid(outcome ~ ., scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Width"))
  )

  print(result_ALL %>%
    group_by(model, outcome) %>%
    dplyr::summarise(Coverage_rate = mean(coverage)) %>%
    ggplot() + geom_bar(aes(x = model, y = Coverage_rate), stat="identity") + facet_grid(outcome ~ .) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", CR"))
  )



  print(result_Z1_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = bias)) + geom_hline(yintercept=0, linetype="dashed", color = "red") + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Bias, Subgroup: Z1 = 1"))
  )


  print(result_Z1_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = square_error)) + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Square Error, Subgroup: Z1 = 1"))
  )

  print(result_Z1_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = AW)) + facet_grid(outcome ~ ., scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Width, Subgroup: Z1 = 1"))
  )

  print(result_Z1_1 %>%
          group_by(model, outcome) %>%
          dplyr::summarise(Coverage_rate = mean(coverage)) %>%
          ggplot() + geom_bar(aes(x = model, y = Coverage_rate), stat="identity") + facet_grid(outcome ~ .) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", CR, Subgroup: Z1 = 1"))
  )


  print(result_Z1_1_Z2_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = bias))+ geom_hline(yintercept=0, linetype="dashed", color = "red") + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Bias, Subgroup: Z1 = 1 & Z2 = 1"))
  )


  print(result_Z1_1_Z2_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = square_error)) + facet_grid(outcome ~ .,  scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Square Error, Subgroup: Z1 = 1 & Z2 = 1"))
  )

  print(result_Z1_1_Z2_1 %>%
          ggplot() + geom_boxplot(aes(x = model, y = AW)) + facet_grid(outcome ~ ., scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", Width, Subgroup: Z1 = 1 & Z2 = 1"))
  )

  print(result_Z1_1_Z2_1 %>%
          group_by(model, outcome) %>%
          dplyr::summarise(Coverage_rate = mean(coverage)) %>%
          ggplot() + geom_bar(aes(x = model, y = Coverage_rate), stat="identity") + facet_grid(outcome ~ .) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Setting ", setting,", CR, Subgroup: Z1 = 1 & Z2 = 1"))
  )

 list(result_ALL, result_Z1_1, result_Z1_1_Z2_1, result_Z1_1_Z2_1_Z3_1)
}



a4 = readRDS("setting1.rds")
a5 = readRDS("setting2.rds")
a6 = readRDS("setting3.rds")


a10_2 = readRDS("setting1_3.rds")
a11_2 = readRDS("setting2_3.rds")
a12_2 = readRDS("setting3_3.rds")

saveRDS(a1, "setting1_final.rds")
saveRDS(a2, "setting2_final.rds")
saveRDS(a12, "setting3_3.rds")

saveRDS(a10, "setting1_mean.rds")
saveRDS(a11, "setting2_mean.rds")
saveRDS(a12, "setting3_mean.rds")


setting = 1
result1 = summary_results(a1, 1)
result1_2 = summary_results(a10_2, 1)

result2 = summary_results(a11, 2)
result2_2 = summary_results(a11_2, 2)


result3 = summary_results(a12, 3)
result3_2 = summary_results(a12_2, 3)


  #


result1[[4]] %>%
  mutate(outcome = factor(outcome, levels = c("continuous", "binary"))) %>%
  group_by(model, outcome) %>%
  dplyr::summarize(absolute_bias = round(abs(mean(bias)), 4) * 100,
                   RMSE = round(sqrt(mean(square_error)), 4)  * 100,
                   AW = round(mean(AW), 4)  * 100,
                   CR = round(mean(coverage), 4) * 100
  ) %>%
  filter(outcome == "continuous")

result1[[3]] %>%
  mutate(outcome = factor(outcome, levels = c("continuous", "binary"))) %>%
  group_by(model, outcome) %>%
  dplyr::summarize(absolute_bias = round(abs(mean(bias)), 4) * 100,
                   RMSE = round(sqrt(mean(square_error)), 4)  * 100,
                   AW = round(mean(AW), 4)  * 100,
                   CR = round(mean(coverage), 4) * 100
  ) %>%
  filter(outcome == "binary")
