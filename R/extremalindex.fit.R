#' Fit a non-stationary regression model to the GPD rate parameter
#' 
#' @name extremalindex.fit
#' 
#' @description This function fits an exponential decay model to empirical estimates of the extremal index of skew surge observations.
#' 
#' @param data A data frame of skew surge observations (variable named \code{skews}).
#' @param run.length This is the number of consecutive non-exceedances between two extreme observations, where we would say they belong to different clusters (as for the standard runs estimate). The default is 10 (i.e., 5 days) but can be inferred from autocorrelation (acf) plots, i.e., the maximum lag before the acf remains close to zero.
#' @param thresh.quantile The quantile of skew surges above which the exponential decay model is required, i.e., the empirical estimates become noisy above this value. This is a single value between 0 and 1. The default is the 0.99-quantile.
#' 
#' @return 2 parameters for the extremal index model fit.
#' 
#' @details \loadmathjax{} The extremalindex.fit consists of two stages. 
#' For values below the threshold (defined by \code{thresh.quantile}), the empirical runs estimate is used as this is smooth over skew surges in this range. 
#' For computational efficiency purposes, these empirical estimates are evaluated on a regular grid of 100 values from the minimum skew surge up to this threshold. 
#' The empirical runs estimate is found using the \code{evd} package.
#' Linear interpolation is used to values between those on the regular grid; the \code{linear.interp} function does this step. 
#' 
#' For skew surges above the threshold, we use the parametric form \mjeqn{\hat\theta(y,r)=\theta-[\theta-\tilde{\theta}(v,r)]\exp\big(-\frac{y-v}{\psi}\big)}{see equation (4.16) in the manuscript} where \mjeqn{y}{y} is the skew surge observation, \mjeqn{v}{v} the threshold, \mjeqn{r}{r} the pre-specified run length, \mjeqn{\tilde{\theta}(v,r)}{tilde(theta)(v,r)} is the empirical runs estimate at the threshold \mjeqn{v}{v} with run length \mjeqn{r}{r}, and \mjeqn{\psi>0, \tilde{\theta}(v,r)\leq\theta\leq1}{psi>0, tilde(theta) <= theta <= 1} are parameters to be estimated.
#' 
#' The parameters are estimated using a weighted-least squares approach with weight \mjeqn{w(y)=\sqrt{c(y,r)-1}}{w(y) = (c(y-r)-1)^(1/2)} where \mjeqn{c(y,r)}{c(y,r)} is the number of clusters above \mjeqn{y}{y} separated by run length \mjeqn{r}{r}.
#' We fit this model using the EI_decaymodel function and optimising via optim.
#' 
#' Further details can be found in Section 4 of \insertCite{DArcy2023;textual}{ESLestimation}. 
#' 
#' @rdname extremalindex.fit
#' 
#' @references \insertRef{DArcy2023}{ESLestimation}
#' 
#' @examples
#' 
#' library(ESLestimation)
#' 
#' # Load in data
#' # This example uses Lowestoft data available within the package 
#' 
#' data(Lowestoft)
#' 
#' # Check the data has the required variables that are correctly named (skews, month, day, maxTide)
#' head(Lowestoft)
#' 
#' # Infer the run length from the acf plot
#' acf(Lowestoft$skews)
#' # Run length of 10 is reasonable
#' 
#' # Check when empirical estimates become noisy
#' # Note this takes a while to run but is an important check
#' library(extRemes)
#' plot(Lowestoft$skews, sapply(Lowestoft$skews, function(y) extremalindex(Lowestoft$skews, y, method='runs', run.length=10)[1]), xlab='Skew surge (m)', ylab=expression(tilde(theta)))
#' abline(v=quantile(Lowestoft$skews, 0.99, names=F))
#' # 0.99-quantile is reasonable
#' 
#' # Fit the extremal index model to skew surges
#' # In this example we use the standard arguments for run.length and thresh.quantile
#' ( rate_par <- rateparam.fit(Lowestoft) )
#' 
#' @export

extremalindex.fit <- function(data, run.length=10, thresh.quantile = 0.99){
  quan <- quantile(data$skews, thresh.quantile, names=F)
  
  EI_lin.interp <- linear.interp(data, thresh.quantile, run.length)
  
  EI_decaymodel <- function(y,par){
    theta <- par[1]
    phi_2 <- par[2]
    
    w <- S <- c()
    if( theta<=1 & theta>=EI_lin.interp(quan) & phi_2>0){
      for(i in 1:length(y)){
        w[i] <- sqrt(extremalindex(data$skews,y[i],method='runs',run.length = run.length)[2]-1)
        S[i] <- w[i]*((1/(extremalindex(data$skews,y[i],method='runs',run.length = run.length)[1]))-(1/(theta-((theta-EI_lin.interp(quan))*exp(-(y[i]-quan)/phi_2)))))^2
      }
      return(sum(S,na.rm=T))
    }
    else{
      return(Inf)
    }
  }
  
  EI_modelfit <- optim(par=c(1,0.5),EI_decaymodel,y=data$skews[which(data$skews>quan)])
  return(EI_modelfit$par)
} 
