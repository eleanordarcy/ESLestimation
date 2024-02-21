monmax.sealevel.dist <- function(z, data, month, year, gpd_par, rate_par, extremalindex_par, q, lin.int, emp){
  
  ## Extract month and year specific tide samples
  tides <- data$maxTide[which(data$month== month & data$year==year)]
  day <- data$day[which(data$month == month & data$year == year)]
  stTides <- data$stTide[which(data$month == month & data$year == year)]
  
  ## Derive sea level by inputting tide series into skew surge distribution function
  F_ss_all <- c()
  if( length(tides)==0 ){
    F_sl <- 1
  }else{
    for(i in 1:length(tides)){
      F_ss_all[i] <- skewsurge.dist(z-tides[i], day[i], month, tides[i], stTides[i], par=c(gpd_par,rate_par), data, q=q, emp = emp)^extremalindex.modelform(z-tides[i], par=extremalindex_par, data, lin.int = lin.int)
    }
    F_sl <- prod(F_ss_all,na.rm=T)
  }
  return(F_sl)
} 
