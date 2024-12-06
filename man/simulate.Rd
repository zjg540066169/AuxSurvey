% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxSurvey.R
\name{simulate}
\alias{simulate}
\title{Simulate Survey Data with Discretized Auxiliary Variables}
\usage{
simulate(N = 3000, discretize = c(3, 5, 10), setting = c(1, 2, 3), seed = NULL)
}
\arguments{
\item{N}{Number of population units to simulate. Default is 3000.}

\item{discretize}{A scale specifying the number of categories for discretizing continuous variables.
The function discretizes both \code{X} and \code{W} into the specified categories. Default is a number among \verb{(3, 5, 10)}.}

\item{setting}{A numeric value to specify the simulation setting. The settings define different relationships
between the outcome variables and the covariates. Possible values are 1, 2, 3, and 4.
Default is a number among \code{c(1, 2, 3)}.}

\item{seed}{An optional random seed for reproducibility. Default is \code{NULL}.}
}
\value{
A list containing two elements:
\itemize{
\item \code{population}: A tibble with the simulated population data, including both continuous and binary outcomes,
as well as auxiliary variables (both raw and discretized).
\item \code{samples}: A tibble with the simulated sample data, where individuals are included based on their estimated
propensity scores.
}
}
\description{
This function simulates survey data with discretized auxiliary variables. It generates a population
dataset with continuous and binary outcomes, and includes auxiliary variables that are discretized
into multiple categories. The function also generates a subset of the population as a sample, based on
the propensity scores.
}
\details{
The function supports multiple simulation settings, where each setting modifies the relationships
between the outcome variables and the covariates.
}
\examples{
# Simulate survey data with setting 1 and discretizing variables 3 categories
data = simulate(N = 3000, discretize = 3, setting = 1, seed = 123)

# Extract population and sample datasets
population = data$population
samples = data$samples

# Examine the simulated population data
head(population)

}
