linear.interp <- function(data, thresh.quantile, run.length){
  Y <- seq(min(data$skews), quantile(data$skews, thresh.quantile, names=F), length=100)
  t <- sapply(Y, function(y) extremalindex(data$skews, y, method = 'runs', run.length = run.length)[1])
  return(approxfun(Y,t))
}
