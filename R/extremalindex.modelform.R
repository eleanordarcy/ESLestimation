extremalindex.modelform <- function(y, par, data, thresh.quantile = 0.99, run.length=10, lin.int){
  ## y is skew surge level
  ## par is parameter output from EI_modelfit
  quan <- quantile(data$skews, thresh.quantile, names=F)
  
  if(y<=quan){
    ## If less than 0.99 quantile, use empirical look up table
    theta <- lin.int(y)
  }
  else{
    ## Use model above
    theta <- par[1]-((par[1]-lin.int(quan))*exp(-(y-quan)/par[2]))
  }
  return(theta)
}
