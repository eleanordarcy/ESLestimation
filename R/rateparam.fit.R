#' Fit a non-stationary regression model to the GPD rate parameter
#' 
#' @name rateparam.fit
#' 
#' @description The function fits a harmonic regression model to the rate parameter of the generalised Pareto distribution (GPD), for skew surge exceedances of a monthly threshold. The model for the rate parameter is expressed through a logit link function, and depends on peak tide and day in year covariates.
#' 
#' @param data A data frame of skew surge observations (named \code{skews}), along with \code{month}, \code{day} and \code{maxTide} observations.
#' @param q The monthly quantile of skew surges to obtain exceedances, this must be between 0 and 1. The default is 0.95 but this should be investigated using standard threshold selection techniques.
#' @param optim.method The method to be used for optimisation of the log-likelihood function for the GPD. The default is \code{'BFGS'} but see \code{\link[stats]{optim}} for details and other options.
#' @param init.par The initial values for parameters to be optimised over in \code{\link[stats]{optim}}. The default are \code{c(0.2,0,100,0.1,0.1)} but these can be altered. 
#' 
#' @return 5 parameters for the GPD rate parameter regression model fit.
#' 
#' @details \loadmathjax{} This function uses maximum likelihood estimation to fit a harmonic-based regression model to the GPD rate parameter, used for skew surge exceedances of a monthly threshold (determined by the quantile \mjeqn{q}{q}).
#' Non stationarity is introduced to the rate parameter of the GPD using two harmonics to capture seasonal variations and skew surge-peak tide dependence. 
#' To fit the model, we use an indicator function to determine whether an observation exceeds its monthly quantile and assume this binary variable follows a Bernoulli distribution.
#' 
#' The rate parameterisation is given by \mjeqn{g(\lambda_{d,t})=g\big(\frac{1-q}{q}\big)+\tilde{d}\beta\sin\big(\frac{2\pi}{f}(d-\phi)\big)+\tilde{x}\big[\alpha^{(x)}+\beta^{(x)}\sin\big(\frac{2\pi}{f}(d-\phi^{(x)})\big)\big]}{See equation (4.11) in the manuscript} for \mjeqn{\alpha^{(x)},\beta,\beta^{(x)},\phi,\phi^{(x)}}{alpha(x), beta, beta(x), phi, phi(x)} parameters to be estimated, \mjeqn{d}{d} day in year (\mjeqn{\tilde{d}\in[-15,15]}{tilde(d)} standardised day in month by subtracting monthly mean), \mjeqn{x}{x} the peak tide observation (\mjeqn{\tilde{x}}{tilde{x}} standardised peak tide by subtracting mean and diving by variance), \mjeqn{f=365}{f=365} the periodicity of the harmonic and \mjeqn{g()}{g} the logit link function. 
#' The Bernoulli log likelihood is given by the function \code{Bernoulli.ll} within \code{rateparam.fit}. This is optimised using \code{\link[stats]{optim}} to obtain parameter estimates for \mjeqn{\lambda}{lambda}.
#' 
#' Further details can be found in Section 4 of \insertCite{DArcy2023;textual}{ESLestimation}. 
#' 
#' @rdname rateparam.fit
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
#' # Fit the haromnic regression model for the rate parameter to an indicator vector
#' # This is 1 for exceedances of the 0.95 monthly quantile, and 0 otherwise.
#' # In this example we use the standard arguments for q, optim.method and init.par
#' ( rate_par <- rateparam.fit(Lowestoft) )
#' 
#' @export



rateparam.fit <- function(data, q=0.95, optim.method='BFGS', init.par=c(0.2,0.2,0.2,0.2,360)){
  data <- data %>%
    group_by(month) %>%
    mutate(stDay = day - mean(day)) %>%
    ungroup()
  
  data$ind <- ifelse(data$skews > ave(data$skews, data$month, FUN = function(x) quantile(x, q)), 1, 0)
  

  Bernoulli_ll<-function(par,z){
    g <- par[1]
    h <- par[2]
    a <- par[3]
    b <- par[4]
    phi <- par[5]
    
    N <- length(z$skews)
    
    f <- 365
    
    beta_1 <- beta_0 <- c()
    for(i in 1:N){
      H <- g*(sin(((2*pi)/f)*(z$day[i]-h)))
      beta_0[i] <- log((1-q)/q)+H*z$stDay[i]
      beta_1[i] <- a+b*(sin(((2*pi)/f)*(z$day[i]-phi)))
    }
    
    lambda <- c() 
    for(i in 1:N){
      lambda[i] <- exp(beta_0[i]+(beta_1[i]*z$stTide[i]))/(1+exp(beta_0[i]+(beta_1[i]*z$stTide[i])))
    }
    
    l<-c()
    for(i in 1:N){
      l[i] <- (z$ind[i]*log(lambda[i]))+((1-z$ind[i])*(log(1-lambda[i])))
    }
    return(-1*sum(l))
  }
  
  Rate_modelfit <- optim(par=init.par,fn=Bernoulli_ll,z=data,method='BFGS')
  return(Rate_modelfit$par)
}
