#' Estimate sea level return levels from the non-stationary model of D`Arcy et al. (2023).
#' 
#' @name returnlevel.est
#' 
#' @description \loadmathjax{} This function brings together the non-stationary GPD model for skew surges with the known tidal regime, to estimate extreme sea levels for a specified annual exceedance probability \mjeqn{p}{p}.
#' 
#' @param p The annual exceedance probability (1/return period) for the required return level estimate. This can be a single value or a vector of probabilities. This must take values between 0 and 1 
#' @param data A data frame of skew surge observations (named \code{skews}), along with \code{month}, \code{day} and \code{maxTide} observations for covariate information.
#' @param gpd_par Scale and shape parameters for the GPD fit to skew surges (found via \code{\link{GPD.fit}} function). This should be a vector of length 5. 
#' @param rate_par Rate parameter for the GPD fit to skew surges (found via \code{\link{rateparam.fit}} function). This should be a vector of length 5.
#' @param extremalindex_par Parameters for the extremal index model fit (found via \code{\link{extremalindex.fit}} function). This should be a vector of length 2.
#' @param gpd.quantile The quantile used to define exceedances and fit the GPD model for skew surges. This is a single value between 0 and 1. The default is 0.95. 
#' @param extremalindex.quantile The quantile of skew surges for the extremal index model, above which the exponential decay model is required, i.e., the empirical estimates become noisy above this value. This is a single value between 0 and 1. The default is the 0.99-quantile.
#' @param run.length This is the number of consecutive non-exceedances between two extreme observations, where we would say they belong to different clusters (as for the standard runs estimate). The default is 10 (i.e., 5 days) but can be inferred from autocorrelation (acf) plots, i.e., the maximum lag before the acf remains close to zero.
#' 
#' @return Return level sea level estimate (in metres).
#' 
#' @details \loadmathjax{} To estimate sea level return levels, we use the derived annual maxima distribution for sea levels of \insertCite{DArcy2023;textual}{ESLestimation}.
#' This uses the fact that sea levels \mjeqn{Z}{Z} can be decomposed into skew surge \mjeqn{Y}{Y} and peak tide \mjeqn{X}{X}. 
#' We detail the methodology for estimating return levels below, with simplified notation but refer the reader to \insertCite{DArcy2023;textual}{ESLestimation} for details.
#' 
#' Skew surges are modelled using non-stationary GPD, for exceedances of a monthly threshold, and a non-stationary monthly empirical distribution for non-exceedances.
#' We account for non-stationarity in terms of seasonality and skew surge-peak tide dependence, using daily, monthly and peak covariates within the model.
#' This model is fit using the \code{\link{GPD.fit}} function (see \code{\link{GPD.fit}} documentation for details) and the output of that function (i.e., the GPD parameter estimates) are used as an input to this function.
#' The rate parameter is separately modelled, with the same covariate information, through the \code{\link{rateparam.fit}} function and its output is also an input to this function.
#' Lastly, we model the temporal dependence of skew surges through a model for the extremal index (see \code{\link{extremalindex.fit}} documentation for details), the output to this function is used as a input to this function.
#' The distribution function is given by \code{skewsurge.dist}.
#' 
#' As tides are deterministic, we chose tidal samples to input into our skew surge distribution to derive a model for the skew surges.
#' This is done via a joint probability methodology, assumed that skew surge-peak tide dependence has been fully captured through the non-stationary skew surge model.
#' So that, \mjeqn{Pr(Z\leq z)=\Pr(X+Y\leq z)=\Pr(Y\leq z-X)=F_{Y}(z-X)}{see equation (4.2) in the manuscript}.
#' Then the annual maxima distribution for sea levels is found by taking the product over the skew surge distribution, for observations and their assoiated covariates within a year.
#' Since skew surges are not independent, use the extremal index model \mjeqn{\hat\theta}{hat(theta)} as follows: \mjeqn{\Pr(M\leq z) = \prod F_Y(z-X)^{\hat\theta}}{see equation (4.2) in the manuscript} for \mjeqn{M}{M} denoting the annual maximum sea level.
#' To account for seasonality, we first derive the monthly maxima distribution using month-specific tidal series.
#' This is given by the \code{monmax.dist} function
#' To account for interannual tidal variations, we do this for each year in the record, using the associated tidal sequence, and then average over all years.
#' The annual maxima distribution function is given by \code{annmax.dist} below.
#' 
#' We use a numerical solver (\code{\link[stats]{uniroot}}) to find the return level \mjeqn{z_p}{z_p} for the associate annual exceedance probablity \mjeqn{p}{p}.
#' That is, to solve the equation \mjeqn{\Pr(M \leq z_p) = 1 - p}{ P(Z < z_p) = 1-p}.
#' Uniroot requires an interval to test over for values of \mjeqn{z_p}{z_p} which are set as the 0.5 quantile of sea levels and the maximum sea level plus 3m for the lower and upper endpoints, respectively.
#' If these endpoints are not wide enough, they will be increased by 0.1m. 
#' 
#' @rdname returnlevel.est
#' 
#' @references \insertRef{DArcy2023}{ESLestimation}
#' 
#' @examples
#' 
#' library(ESLestimation)
#' 
#' # Load in data
#' # This example uses Lowestoft data available within the package 
#' data(Lowestoft)
#' 
#' # Check the data has the required variables that are correctly named (skews, month, day, maxTide)
#' head(Lowestoft)
#' 
#' # Estimate 10 year return level (annual exceedance probability p=0.1)
#' # Use model fit functions to get parameter estimates for the GPD, rate and extremal index.
#' # Default quantiles and run length are used.
#' # This takes ~18 minutes to run for a single probability p.
#' ( RL <- returnlevel.est(0.1, Lowestoft, gpd_par=GPD.fit(Lowestoft), rate_par=rateparam.fit(Lowestoft), extremalindex_par=extremalindex.fit(Lowestoft)) )
#' 
#' @export
returnlevel.est <- function(p, data, gpd_par, rate_par, extremalindex_par, gpd.quantile = 0.95, extremalindex.quantile = 0.99, run.length = 10){
  period <- (max(data$year)+1)-min(data$year)
  base_yr <- min(data$year)
  return_levels <- c()
  
  EI_lin.interp <- linear.interp(data, extremalindex.quantile, run.length)
  empiricals <- get.empiricals(data)
  
  
  p <- c(p) 
  for(j in 1:length(p)){
    print(paste0('Estimating ', 1/p[j], ' year return level'))
    annmax.sealevel.dist <- function(z){
      F_k <- c()
      for(i in 1:period){
        F_m_s <- c()
        F_m_s[1] <- monmax.sealevel.dist(z,data,1,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[2] <- monmax.sealevel.dist(z,data,2,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[3] <- monmax.sealevel.dist(z,data,3,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[4] <- monmax.sealevel.dist(z,data,4,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[5] <- monmax.sealevel.dist(z,data,5,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[6] <- monmax.sealevel.dist(z,data,6,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[7] <- monmax.sealevel.dist(z,data,7,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[8] <- monmax.sealevel.dist(z,data,8,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[9] <- monmax.sealevel.dist(z,data,9,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[10] <- monmax.sealevel.dist(z,data,10,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[11] <- monmax.sealevel.dist(z,data,11,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        F_m_s[12] <- monmax.sealevel.dist(z,data,12,base_yr + (i-1), gpd_par, rate_par, extremalindex_par, q=gpd.quantile, lin.int = EI_lin.interp, emp = empiricals)
        
        F_k[i] <- prod(F_m_s)
      }
      mean(F_k)-(1-p[j])
    }

    solver.interval.min <- quantile(data$skews+data$maxTide, 0.5, names=F)
    solver.interval.max <- max(data$skews+data$maxTide)+3
    while( sign(annmax.sealevel.dist(solver.interval.min)) == sign(annmax.sealevel.dist(solver.interval.max)) ){
      solver.interval.min <- solver.interval.min - 0.1
      solver.interval.max <- solver.interval.max + 0.1
    }
    return_levels[j] <- uniroot(annmax.sealevel.dist, c(solver.interval.min, solver.interval.max))$root
    
    return(return_levels)
  }
}
