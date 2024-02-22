get.data.unif.margins <- function(data, par, emp, q){
  N <- length(data$skews)
  U <- c()
  for(i in 1:N){
    U[i]<- skewsurge.dist(data$skews[i], data$day[i], data$month[i], data$maxTide[i], data$stTide[i], par=par, data=data, q=q, emp=emp)
  }
  
  return(U)
}