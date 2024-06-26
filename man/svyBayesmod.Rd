% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxSurvey.R
\name{svyBayesmod}
\alias{svyBayesmod}
\title{Title}
\usage{
svyBayesmod(
  svysmpl,
  svypopu,
  outcome_formula,
  BayesFun,
  subset = NULL,
  family = gaussian(),
  invlvls,
  weights = NULL,
  nskip = 1000,
  npost = 1000,
  nchain = 4,
  printmod = TRUE,
  doFigure = FALSE,
  useTrueSample = F,
  stan_verbose = F,
  shortest_CI = F
)
}
\arguments{
\item{svysmpl}{A dataframe or tibble to represent 'samples'.}

\item{svypopu}{A dataframe or tibble to represent 'population'.}

\item{outcome_formula}{A formula for stan.}

\item{BayesFun}{Function name for stan}

\item{subset}{A character vector. Each element is a string representing a filtering condition to select subset of samples and population. Default is NULL. When this parameter is NULL, the analysis is only performed on the whole data. If subsets are specified, the estimates for the whole data will also be calculated.}

\item{family}{The distribution family of outcome variable. Currently we only support \code{\link[stats]{gaussian}} and \code{\link[stats]{binomial}}.}

\item{invlvls}{A numeric vector of values. Each specifies a confidence level of CI for estimators. If more than one values are specified, then multiple CIs are calculated.}

\item{weights}{A numeric vector of case weights. The length should be equal to the number of cases in 'samples'.}

\item{nskip}{A integer to specify the number of burn-in iterations of each chain in MCMC for stan models. Default is 1000.}

\item{npost}{A integer to specify the number of posterior sampling iteration of each chain in MCMC for stan models. Default is 1000.}

\item{nchain}{A integer to specify the number of MCMC chains for stan models. Default is 4.}

\item{printmod}{A logical to indicate if print posterior estimates}

\item{doFigure}{A logical to indicate if print MCMC figures}

\item{useTrueSample}{A logical to indicate if the estimator uses samples information}

\item{stan_verbose}{A logical to indicate if print MCMC information in stan}

\item{shortest_CI}{A logical scalar; if true, the calculated credible intervals for stan models are highest posterior density intervals. Otherwise the intervals are symmetric. Default is false.}
}
\value{
A list. Each element contains the Bayesian estimate and CIs for a subset or the whole data analysis.
}
\description{
Title
}
