% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI.est.R
\name{CI.est}
\alias{CI.est}
\title{Confidence intervals for sea level return levels from the non-stationary model of D`Arcy et al. (2023).}
\usage{
CI.est(
  p,
  data,
  ci.prob,
  gpd_par,
  rate_par,
  extremalindex_par,
  block.length = 10,
  n.boot = 200,
  gpd.quantile = 0.95,
  optim.method = "BFGS",
  EI.run.length = 10,
  EI.quantile = 0.99
)
}
\arguments{
\item{p}{The annual exceedance probability (1/return period) for the required return level estimate. This can be a single value or a vector of probabilities. This must take values between 0 and 1.}

\item{data}{A data frame of skew surge observations (named \code{skews}), along with \code{month}, \code{day} and \code{maxTide} observations for covariate information.}

\item{ci.prob}{The probability for the confidence interval width. Default is 0.95. This should be a value between 0 and 1.}

\item{gpd_par}{Scale and shape parameters for the GPD fit to skew surges (found via \code{\link{GPD.fit}} function). This should be a vector of length 5.}

\item{rate_par}{Rate parameter for the GPD fit to skew surges (found via \code{\link{rateparam.fit}} function). This should be a vector of length 5.}

\item{extremalindex_par}{Parameters for the extremal index model fit (found via \code{\link{extremalindex.fit}} function). This should be a vector of length 2.}

\item{block.length}{The block length for the stationary bootstrap procedure, and should represent the approximate duration of a storm. The default is 10, corresponding to approximately 5 days.}

\item{n.boot}{The number of bootstrap samples to use to obtain confidence intervals. The default is 200.}

\item{gpd.quantile}{The quantile used to define exceedances and fit the GPD model for skew surges. This is a single value between 0 and 1. The default is 0.95.}

\item{optim.method}{The method for optimisation when refitting the models to bootstrap samples. The default is \code{'BFGS'}, see \code{\link[statas]{optim}} for more details.}

\item{EI.run.length}{This is the number of consecutive non-exceedances between two extreme observations, where we would say they belong to different clusters (as for the standard runs estimate). This is used for fitting the extremal index model. The default is 10 (i.e., 5 days) but can be inferred from autocorrelation (acf) plots, i.e., the maximum lag before the acf remains close to zero.}

\item{EI.quantile}{The quantile of skew surges for the extremal index model, above which the exponential decay model is required, i.e., the empirical estimates become noisy above this value. This is a single value between 0 and 1. The default is the 0.99-quantile.}
}
\value{
(1-\code{ci.prob}/2)), 0.5, \code{ci.prob}+(1-\code{ci.prob}/2)) quantiles of return level sea level estimates (in metres) over the specified number of bootstrap samples.
}
\description{
This function gives confidence intervals on the sea level return level estimates from the non-stationary model of \insertCite{DArcy2023;textual}{ESLestimation}.
}
\details{
\loadmathjax{} A stationary bootstrap procedure is used to obtain confidence intervals on the sea level return level estimates given by the \code{\link{returnlevel.est}} function.
The block length is simulated from a Geometric distribution with mean as the reciprocal of the input \code{block.length}. 
This can be inferred from an autocorrelation function plot.

We account for uncertainty at each stage of the modelling procedure by recalculating thresholds, re-estimating model parameters and the empirical distribution for each bootstrap sample.
Once the required number of bootstrap samples is reached (\code{n.boot}), we find the required quantiles (indicated by \code{ci.prob}) for confidence intervals.
The 0.5 quantile is also returned.
}
\examples{

library(ESLestimation)

# Load in data
# This example uses Lowestoft data available within the package 
data(Lowestoft)

# Check the data has the required variables that are correctly named (skews, month, day, maxTide)
head(Lowestoft)

# Estimate 10 year return level (annual exceedance probability p=0.1) and it's associated 95\\% confidence interval.
# Use model fit functions to get parameter estimates for the GPD, rate and extremal index.
# Default parameters are used everywhere - this gives a 95\% confidence interval
CI.est(0.1, Lowestoft, ci.prob, gpd_par = GPD.fit(Lowestoft), rate_par = rateparam.fit(Lowestoft), extremalindex_par = extremalindex.fit(Lowestoft), block.length=10, n.boot=200, gpd.quantile=0.95, optim.method='BFGS', EI.run.length=10, EI.quantile = 0.99)

}
\references{
\insertRef{DArcy2023}{ESLestimation}
}
