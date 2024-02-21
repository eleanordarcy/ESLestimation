skewsurge.dist <- function(y, day, month, tide, stTide, par, data, q, emp){
  ## par is the parameter vector for the scale, shape and rate parameters (a,b,p,c,b1,p1,xi)
  ## a,b,p characterise the harmonic on the scale parameter
  ## c is the linear trend for tide on the scale parameter
  ## 5 rate parameters characterise the harmonic in model R1 fit above
  ## xi is the shape parameter
  ## q is the quantile for the GPD threshold
  
  a <- par[1]
  b <- par[2]
  p <- par[3]
  c <- par[4]
  
  rate <- par[6:10]
  
  xi <- par[5]
  
  i <- day
  j <- month
  t <- tide
  st <- stTide
  
  if(y <= quantile(data$skews[which(data$month==j)], q, names=F)){
    return(emp.dist(data, y, obs.tide=t, mnth=j, emp = emp))
  }
  else{
    ## Scale harmonic
    sigma <- a+b*(sin(((2*pi)/365)*(i-p)))+t*c
    
    ## Rate harmonic
    sd <- i-mean(data$day[which(data$month==j)])
    beta_0 <- log((1-q)/q)+sd*(rate[1]*sin(((2*pi)/365)*(i-rate[2])))
    beta_1 <- rate[3]+rate[4]*(sin(((2*pi)/365)*(i-rate[5])))
    exprob <- exp(beta_0+(beta_1*st))/(1+exp(beta_0+(beta_1*st)))
    
    return( 1-(exprob*(1+xi*((y-quantile(data$skews[which(data$month==j)], q, names=F))/sigma))^(-1/xi)))
  }
}
