#' Fit the non-stationary GPD model to skew surges
#' 
#' @name GPD.fit
#' 
#' @description The function fits a non-stationary generalised Pareto distribution (GPD) to skew surge exceedances of a monthly threshold. This has a non-stationary scale parameter to represent seasonal variability (via a harmonic) and skew surge-peak tide dependence (via a linear trend).
#' 
#' @param data A data frame of skew surge observations (named \code{skews}), along with \code{month}, \code{day} and \code{maxTide} observations.
#' @param q The monthly quantile of skew surges to obtain exceedances and must be between 0 and 1. The default is 0.95 but this should be investigated using standard threshold selection techniques.
#' @param optim.method The method to be used for optimisation of the log-likelihood function for the GPD. The default is \code{'BFGS'} but see \code{\link[stats]{optim}} for details and other options.
#' @param init.par The initial values for parameters to be optimised over in \code{\link[stats]{optim}}. The default are \code{c(0.2,0,100,0.1,0.1)} but these can be altered. 
#' 
#' @return Parameters for the GPD model fit: 4 parameters characterising the scale parameter and 1 shape parameter.
#' 
#' @details \loadmathjax{} This function uses maximum likelihood estimation to fit a non-stationary GPD to skew surge exceedances of a monthly threshold (determined by the quantile \mjeqn{q}{q}).
#' Non stationarity is introduced to the scale parameter of the GPD using a harmonic to capture seasonal variations and a linear trend for tide, to capture skew surge-peak tide dependence. 
#' This parametrisation is given by \mjeqn{\sigma_{d,t}=\alpha+\beta\sin[(2\pi/f)(d-\phi)]+\gamma x}{sigma_{d,t}=alpha + beta*sin[((2pi)/f)*(d-phi)]+gamma*x} for \mjeqn{\alpha,\beta,\phi,\gamma}{alpha, beta, phi, gamma} parameters to be estimated, \mjeqn{x}{x} the peak tide observation, \mjeqn{d\in[1:365]}{d} the day in year covariate and \mjeqn{f=365}{f=365} the periodicity of the harmonic.
#' 
#' The GPD log likelihood is given by the function \code{GPD_ll} within \code{GPD.fit} and then optimised using \code{\link[stats]{optim}} to obtain parameter estimates for \mjeqn{\sigma}{sigma} (i.e., \mjeqn{\alpha,\beta,\phi,\gamma}{alpha, beta, phi, gamma}) and \mjeqn{\xi}{xi} the GPD shape parameter.
#' 
#' A Normal prior distribution is incorporated onto the log likelihood to reduce uncertainty associated with shape parameter estimation. 
#' This prior information is based on spatial information about shape parameter estimates across the UK - see the Coastal Flood Boundary report (2018).
#' 
#' Further details can be found in Section 4 of \insertCite{DArcy2023;textual}{ESLestimation}.

#' @rdname GPD.fit
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
#' # Fit the non-stationary GPD model to exceedances of the 0.95 monthly quantile.
#' # In this example we use the standard arguments for q, optim.method and init.par
#' ( gpd_par <- GPD.fit(Lowestoft) )
#' 
#' @export

GPD.fit <- function(data, q=0.95, optim.method='BFGS', init.par=c(0.2,0,100,0.1,0.1)){
  data_ex <- data %>%
    group_by(month) %>%
    mutate(skews = skews - quantile(skews, q)) %>%
    filter(skews >= 0) %>%
    ungroup()
  
  GPD_ll<-function(par,z){
    
    a <- par[1]
    b <- par[2]
    p <- par[3]
    d <- par[4]
    xi <- par[5]
    
    f <- 365 # 1 year
    
    N <- length(z$skews)
    
    sig<-c()
    for(i in 1:N){
      sig[i]<-a+b*(sin(((2*pi)/f)*(z$day[i]-p)))+d*z$maxTide[i]
    }
    
    if( (all(sig>0)) ){
      
      T<-c()
      for(i in 1:N){
        if( all(1+xi*((z$skews[i])/sig[i])>0) ){
          T[i]<-TRUE
        }else{
          T[i]<-FALSE
        }
      }
      
      if( all(T==TRUE) ){
        l<-c()
        for(i in 1:N){
          l[i]<-log(1+xi*(z$skews[i]/sig[i]))
        }
        return(-1*sum(log(sig))-((1/xi+1)*(sum(l)))+dnorm( x=xi, mean=0.0119, sd = 0.03432, log=T))
      }
      
      else if( round(xi,3)==0 ){
        k<-c()
        for(i in 1:N){
          k[i]<-sum(z$skews[i]/sig[i])
        }
        return(-1*log(sig)-sum(k)+dnorm( x=xi, mean=0.0119, sd = 0.03432, log=T))
      }
      
      else{
        return(-1e7)
      }
    }
    else{
      return(-1e6)
    }
  }
  
  GPD_modelfit <- optim(par=init.par, GPD_ll, z=data_ex, control=list(fnscale=-1), method=optim.method)
  
  return(GPD_modelfit$par)
}
