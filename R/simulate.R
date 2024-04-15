library(tidyverse)
library(mgcv)

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
    pi = invlogit(- 1.25 - Z[,1] + 1.25*Z[,2] - 0.75*Z[,3] + 0.75*X - 0.10*X^2)
    Y1 = rnorm(N, 15 + 2.5*Z[,1] - Z[,2] + Z[,3] - 2*X + 3.75*X^2, 3)
    Y2 = rbinom(N, 1, invlogit(-2.5 + 0.75*Z[,1] - 2.5*Z[,2] + 1.5*Z[,3] - 0.25*X + 1.5*X^2))
  }
  if(setting == 2){
    pi = invlogit(- 1.25 - Z[,1] + 1.25*Z[,2] - 0.75*Z[,3] + 0.75*W - 0.10*W^2)
    Y1 = rnorm(N, 15 + 2.5*Z[,1] - Z[,2] + Z[,3] - 2*X + 3.75*X^2, 3)
    Y2 = rbinom(N, 1, invlogit(-2.5 + 0.75*Z[,1] - 2.5*Z[,2] + 1.5*Z[,3] - 0.25*X + 1.5*X^2))
  }
  if(setting == 3){
    pi = invlogit(-0.9 - 0.5 * Z[,1] + 0.75 * Z[,2] - Z[,3] + 0.5 * X - 0.05 * X^2 + 0.5 * Z[,1] * X - 0.75 * Z[,1] * X^2)
    Y1 = rnorm(N, 15 + 2.5 * Z[,1] - Z[,2] + Z[,3] - 2 * X + X^2 + Z[,1] * X - 2.5 * Z[,1] * X^2, 2)
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
    logit_true_pi = gtools::logit(pi)
  )

  U = runif(N)
  population$inclusion = ifelse(pi > U, T, F)
  samples = population[population$inclusion == T,]
  colnames(population) = c("id", "Z1", "Z2", "Z3", "X", "W", colnames(auX), colnames(auW), "Y1", "Y2", "true_pi", "logit_true_pi", "inclusion")
  colnames(samples) = c("id", "Z1", "Z2", "Z3", "X", "W", colnames(auX), colnames(auW), "Y1", "Y2", "true_pi", "logit_true_pi", "inclusion")

  ps_model = BART3::gbart(as.matrix(population[, c("Z1", "Z2", "Z3", "X")]), population$inclusion, type = "pbart", ntree=50, keepevery= 1, ndpost = 100)
  population$estimated_pi = predict(ps_model, as.matrix(population[, c("Z1", "Z2", "Z3", "X")]))$prob.test.mean
  samples$estimated_pi = predict(ps_model, as.matrix(samples[, c("Z1", "Z2", "Z3", "X")]))$prob.test.mean
  population$logit_estimated_pi = logit(population$estimated_pi)
  samples$logit_estimated_pi = logit(samples$estimated_pi)

  ps_model_with_W <- BART3::gbart(as.matrix(population[, c("Z1", "Z2", "Z3", "X", "W")]), population$inclusion, type = "pbart", ntree=50, keepevery= 1, ndpost = 100)
  population$estimated_pi_with_W = predict(ps_model_with_W, as.matrix(population[, c("Z1", "Z2", "Z3", "X", "W")]))$prob.test.mean
  samples$estimated_pi_with_W = predict(ps_model_with_W, as.matrix(samples[, c("Z1", "Z2", "Z3", "X", "W")]))$prob.test.mean
  population$logit_estimated_pi_with_W = logit(population$estimated_pi_with_W)
  samples$logit_estimated_pi_with_W = logit(samples$estimated_pi_with_W)

  samples$inclusion = NULL
  population$inclusion = NULL
  return(list(population = as_tibble(population), samples = as_tibble(samples)))
}
