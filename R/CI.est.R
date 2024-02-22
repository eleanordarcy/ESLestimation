#' Confidence intervals for sea level return levels from the non-stationary model of D`Arcy et al. (2023).
#' 
#' @name CI.est
#' 
#' @description This function gives confidence intervals on the sea level return level estimates from the non-stationary model of \insertCite{DArcy2023;textual}{ESLestimation}.
#' 
#' @param p The annual exceedance probability (1/return period) for the required return level estimate. This can be a single value or a vector of probabilities. This must take values between 0 and 1. 
#' @param data A data frame of skew surge observations (named \code{skews}), along with \code{month}, \code{day} and \code{maxTide} observations for covariate information.
#' @param ci.prob The probability for the confidence interval width. Default is 0.95. This should be a value between 0 and 1.
#' @param gpd_par Scale and shape parameters for the GPD fit to skew surges (found via \code{\link{GPD.fit}} function). This should be a vector of length 5. 
#' @param rate_par Rate parameter for the GPD fit to skew surges (found via \code{\link{rateparam.fit}} function). This should be a vector of length 5.
#' @param extremalindex_par Parameters for the extremal index model fit (found via \code{\link{extremalindex.fit}} function). This should be a vector of length 2.
#' @param block.length The block length for the stationary bootstrap procedure, and should represent the approximate duration of a storm. The default is 10, corresponding to approximately 5 days.
#' @param n.boot The number of bootstrap samples to use to obtain confidence intervals. The default is 200.
#' @param gpd.quantile The quantile used to define exceedances and fit the GPD model for skew surges. This is a single value between 0 and 1. The default is 0.95. 
#' @param optim.method The method for optimisation when refitting the models to bootstrap samples. The default is \code{'BFGS'}, see \code{\link[statas]{optim}} for more details.
#' @param EI.run.length This is the number of consecutive non-exceedances between two extreme observations, where we would say they belong to different clusters (as for the standard runs estimate). This is used for fitting the extremal index model. The default is 10 (i.e., 5 days) but can be inferred from autocorrelation (acf) plots, i.e., the maximum lag before the acf remains close to zero.
#' @param EI.quantile The quantile of skew surges for the extremal index model, above which the exponential decay model is required, i.e., the empirical estimates become noisy above this value. This is a single value between 0 and 1. The default is the 0.99-quantile.
#' 
#' @return (1-\code{ci.prob}/2)), 0.5, \code{ci.prob}+(1-\code{ci.prob}/2)) quantiles of return level sea level estimates (in metres) over the specified number of bootstrap samples.
#' 
#' @details \loadmathjax{} A stationary bootstrap procedure is used to obtain confidence intervals on the sea level return level estimates given by the \code{\link{returnlevel.est}} function.
#' The block length is simulated from a Geometric distribution with mean as the reciprocal of the input \code{block.length}. 
#' This can be inferred from an autocorrelation function plot.
#' 
#' We account for uncertainty at each stage of the modelling procedure by recalculating thresholds, re-estimating model parameters and the empirical distribution for each bootstrap sample.
#' Once the required number of bootstrap samples is reached (\code{n.boot}), we find the required quantiles (indicated by \code{ci.prob}) for confidence intervals.
#' The 0.5 quantile is also returned.
#' 
#' @rdname CI.est
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
#' # Estimate 10 year return level (annual exceedance probability p=0.1) and it's associated 95\% confidence interval.
#' # Use model fit functions to get parameter estimates for the GPD, rate and extremal index.
#' # Default parameters are used everywhere - this gives a 95% confidence interval
#' CI.est(0.1, Lowestoft, ci.prob, gpd_par = GPD.fit(Lowestoft), rate_par = rateparam.fit(Lowestoft), extremalindex_par = extremalindex.fit(Lowestoft), block.length=10, n.boot=200, gpd.quantile=0.95, optim.method='BFGS', EI.run.length=10, EI.quantile = 0.99)
#' 
#' @export

CI.est <- function(p, data, ci.prob, gpd_par, rate_par, extremalindex_par, block.length=10, n.boot=200, gpd.quantile=0.95, optim.method='BFGS',
                   EI.run.length=10, EI.quantile = 0.99){
  
  EI_lin.interp <- linear.interp(data)
  empiricals <- get.empiricals(data)
  

  RL <- c()
  
  N <- length(data$year)
  
  data.unif <- get.data.unif.margins(data=data, par=c(gpd_par, rate_par), emp=empiricals, q=gpd.quantile)
  
  set.seed(111)
  for(j in 1:n.boot){
    print(paste0('Bootstrap: ',j))
    
    boot <- stat.boot(data.unif, block.length)
    
    data <- data %>%
      group_by(month) %>%
      mutate(stDay = day - mean(day)) %>%
      ungroup()
    
    data_boot <- data.frame(data$year, data$month, data$day, data$maxTide, data$stTide, data$stDay, data$date, boot)
    names(data_boot) <- c('year','month','day','maxTide','stTide','stDay','date','unif')
    
    data_boot$skews <- NULL
    for(i in 1:N){
      
      if( data_boot$unif[i] <= gpd.quantile){
        obvs_m <- data %>% dplyr::filter( data$month == data_boot$month[i])
        q1 <- quantile(obvs_m$maxTide, 0.333, names=F)
        q2 <- quantile(obvs_m$maxTide, 0.667, names=F)
        
        t <- data_boot$maxTide[i]
        if(t<=q1){
          obvs_m_q1 <- obvs_m %>% dplyr::filter( maxTide <= q1)
          s <- sort(obvs_m_q1$skews)
          data_boot$skews[i] <- s[ceiling.custom(length(obvs_m_q1$year)*data_boot$unif[i])]
        }
        else if( t>q1 & t<=q2){
          obvs_m_q2 <- obvs_m %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
          s<- sort(obvs_m_q2$skews)
          data_boot$skews[i] <- s[ceiling.custom(length(obvs_m_q2$year)*data_boot$unif[i])]
        }
        else{
          obvs_m_q3 <- obvs_m %>% dplyr::filter( maxTide > q2)
          s<- sort(obvs_m_q3$skews)
          data_boot$skews[i] <- s[ceiling.custom(length(obvs_m_q3$year)*data_boot$unif[i])]
        }
      }
      
      else{
        beta_0 <- log((1-gpd.quantile)/gpd.quantile)+data_boot$stDay[i]*(rate_par[1]*sin(((2*pi)/365)*(data_boot$day[i]-rate_par[2])))
        beta_1 <- rate_par[3]+rate_par[4]*(sin(((2*pi)/365)*(data_boot$day[i]-rate_par[5])))
        exprob <- exp(beta_0+(beta_1*data_boot$stTide[i]))/(1+exp(beta_0+(beta_1*data_boot$stTide[i])))
        
        
        sigma <- gpd_par[1]+gpd_par[2]*sin(((2*pi)/365)*(data_boot$day[i]-gpd_par[3]))+gpd_par[4]*data_boot$maxTide[i]
        xi <- gpd_par[5]
        
        m <- data_boot$month[i]
        data_boot$skews[i] <- quantile(data$skews[which(data$month==m)], gpd.quantile) + (sigma/xi)*(((1-data_boot$unif[i])/exprob)^(-xi)-1)
      }
    }
    plot(data_boot$date, data_boot$skews)
    
    Jan <- data_boot %>% dplyr::filter(month==1)
    Feb <- data_boot %>% dplyr::filter(month==2)
    Mar <- data_boot %>% dplyr::filter(month==3)
    Apr <- data_boot %>% dplyr::filter(month==4)
    May <- data_boot %>% dplyr::filter(month==5)
    Jun <- data_boot %>% dplyr::filter(month==6)
    Jul <- data_boot %>% dplyr::filter(month==7)
    Aug <- data_boot %>% dplyr::filter(month==8)
    Sep <- data_boot %>% dplyr::filter(month==9)
    Oct <- data_boot %>% dplyr::filter(month==10)
    Nov <- data_boot %>% dplyr::filter(month==11)
    Dec <- data_boot %>% dplyr::filter(month==12)

    Z <- seq(0-max(Jan$maxTide),7-min(Jan$maxTide),by=0.001) #range of skews
    F1 <- F2 <- F3 <- c()
    q1 <- sort(Jan$maxTide)[0.333*length(Jan$maxTide)]
    q2 <- sort(Jan$maxTide)[0.667*length(Jan$maxTide)]
    Jan1 <- Jan %>% dplyr::filter( maxTide <= q1)
    Jan2 <- Jan %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Jan3 <- Jan %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i] <- sum(Jan1$skews<=Z[i])/length(Jan1$skews)
      F2[i] <- sum(Jan2$skews<=Z[i])/length(Jan2$skews)
      F3[i] <- sum(Jan3$skews<=Z[i])/length(Jan3$skews)
    }
    Jan1_emp_skews <- lut(outputs=round(F1,3),inputs=round(Z,3))
    Jan2_emp_skews <- lut(outputs=round(F2,3),inputs=round(Z,3))
    Jan3_emp_skews <- lut(outputs=round(F3,3),inputs=round(Z,3))

    Z <- seq(0-max(Feb$maxTide),7-min(Feb$maxTide),by=0.001) #range of skews
    F1 <- F2<-F3<-c()
    q1 <- sort(Feb$maxTide)[0.333*length(Feb$maxTide)]
    q2 <- sort(Feb$maxTide)[0.667*length(Feb$maxTide)]
    Feb1 <- Feb %>% dplyr::filter( maxTide <= q1)
    Feb2 <- Feb %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Feb3 <- Feb %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i] <- sum(Feb1$skews<=Z[i])/length(Feb1$skews)
      F2[i] <- sum(Feb2$skews<=Z[i])/length(Feb2$skews)
      F3[i] <- sum(Feb3$skews<=Z[i])/length(Feb3$skews)
    }
    Feb1_emp_skews <- lut(outputs=round(F1,3),inputs=round(Z,3))
    Feb2_emp_skews <- lut(outputs=round(F2,3),inputs=round(Z,3))
    Feb3_emp_skews <- lut(outputs=round(F3,3),inputs=round(Z,3))


    Z<-seq(0-max(Mar$maxTide),7-min(Mar$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Mar$maxTide)[0.333*length(Mar$maxTide)]
    q2<-sort(Mar$maxTide)[0.667*length(Mar$maxTide)]
    Mar1 <- Mar %>% dplyr::filter( maxTide <= q1)
    Mar2 <- Mar %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Mar3 <- Mar %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Mar1$skews<=Z[i])/length(Mar1$skews)
      F2[i]<-sum(Mar2$skews<=Z[i])/length(Mar2$skews)
      F3[i]<-sum(Mar3$skews<=Z[i])/length(Mar3$skews)
    }
    Mar1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Mar2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Mar3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Apr$maxTide),7-min(Apr$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Apr$maxTide)[0.333*length(Apr$maxTide)]
    q2<-sort(Apr$maxTide)[0.667*length(Apr$maxTide)]
    Apr1 <- Apr %>% dplyr::filter( maxTide <= q1)
    Apr2 <- Apr %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Apr3 <- Apr %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Apr1$skews<=Z[i])/length(Apr1$skews)
      F2[i]<-sum(Apr2$skews<=Z[i])/length(Apr2$skews)
      F3[i]<-sum(Apr3$skews<=Z[i])/length(Apr3$skews)
    }
    Apr1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Apr2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Apr3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(May$maxTide),7-min(May$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(May$maxTide)[0.333*length(May$maxTide)]
    q2<-sort(May$maxTide)[0.667*length(May$maxTide)]
    May1 <- May %>% dplyr::filter( maxTide <= q1)
    May2 <- May %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    May3 <- May %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(May1$skews<=Z[i])/length(May1$skews)
      F2[i]<-sum(May2$skews<=Z[i])/length(May2$skews)
      F3[i]<-sum(May3$skews<=Z[i])/length(May3$skews)
    }
    May1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    May2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    May3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))


    Z<-seq(0-max(Jun$maxTide),7-min(Jun$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Jun$maxTide)[0.333*length(Jun$maxTide)]
    q2<-sort(Jun$maxTide)[0.667*length(Jun$maxTide)]
    Jun1 <- Jun %>% dplyr::filter( maxTide <= q1)
    Jun2 <- Jun %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Jun3 <- Jun %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Jun1$skews<=Z[i])/length(Jun1$skews)
      F2[i]<-sum(Jun2$skews<=Z[i])/length(Jun2$skews)
      F3[i]<-sum(Jun3$skews<=Z[i])/length(Jun3$skews)
    }
    Jun1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Jun2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Jun3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Jul$maxTide),7-min(Jul$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Jul$maxTide)[0.333*length(Jul$maxTide)]
    q2<-sort(Jul$maxTide)[0.667*length(Jul$maxTide)]
    Jul1 <- Jul %>% dplyr::filter( maxTide <= q1)
    Jul2 <- Jul %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Jul3 <- Jul %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Jul1$skews<=Z[i])/length(Jul1$skews)
      F2[i]<-sum(Jul2$skews<=Z[i])/length(Jul2$skews)
      F3[i]<-sum(Jul3$skews<=Z[i])/length(Jul3$skews)
    }
    Jul1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Jul2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Jul3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Aug$maxTide),7-min(Aug$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Aug$maxTide)[0.333*length(Aug$maxTide)]
    q2<-sort(Aug$maxTide)[0.667*length(Aug$maxTide)]
    Aug1 <- Aug %>% dplyr::filter( maxTide <= q1)
    Aug2 <- Aug %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Aug3 <- Aug %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Aug1$skews<=Z[i])/length(Aug1$skews)
      F2[i]<-sum(Aug2$skews<=Z[i])/length(Aug2$skews)
      F3[i]<-sum(Aug3$skews<=Z[i])/length(Aug3$skews)
    }
    Aug1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Aug2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Aug3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Sep$maxTide),7-min(Sep$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Sep$maxTide)[0.333*length(Sep$maxTide)]
    q2<-sort(Sep$maxTide)[0.667*length(Sep$maxTide)]
    Sep1 <- Sep %>% dplyr::filter( maxTide <= q1)
    Sep2 <- Sep %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Sep3 <- Sep %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Sep1$skews<=Z[i])/length(Sep1$skews)
      F2[i]<-sum(Sep2$skews<=Z[i])/length(Sep2$skews)
      F3[i]<-sum(Sep3$skews<=Z[i])/length(Sep3$skews)
    }
    Sep1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Sep2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Sep3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Oct$maxTide),7-min(Oct$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Oct$maxTide)[0.333*length(Oct$maxTide)]
    q2<-sort(Oct$maxTide)[0.667*length(Oct$maxTide)]
    Oct1 <- Oct %>% dplyr::filter( maxTide <= q1)
    Oct2 <- Oct %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Oct3 <- Oct %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Oct1$skews<=Z[i])/length(Oct1$skews)
      F2[i]<-sum(Oct2$skews<=Z[i])/length(Oct2$skews)
      F3[i]<-sum(Oct3$skews<=Z[i])/length(Oct3$skews)
    }
    Oct1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Oct2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Oct3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Nov$maxTide),7-min(Nov$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Nov$maxTide)[0.333*length(Nov$maxTide)]
    q2<-sort(Nov$maxTide)[0.667*length(Nov$maxTide)]
    Nov1 <- Nov %>% dplyr::filter( maxTide <= q1)
    Nov2 <- Nov %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Nov3 <- Nov %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Nov1$skews<=Z[i])/length(Nov1$skews)
      F2[i]<-sum(Nov2$skews<=Z[i])/length(Nov2$skews)
      F3[i]<-sum(Nov3$skews<=Z[i])/length(Nov3$skews)
    }
    Nov1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Nov2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Nov3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    Z<-seq(0-max(Dec$maxTide),7-min(Dec$maxTide),by=0.001) #range of skews
    F1<-F2<-F3<-c()
    q1<-sort(Dec$maxTide)[0.333*length(Dec$maxTide)]
    q2<-sort(Dec$maxTide)[0.667*length(Dec$maxTide)]
    Dec1 <- Dec %>% dplyr::filter( maxTide <= q1)
    Dec2 <- Dec %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
    Dec3 <- Dec %>% dplyr::filter( maxTide >q2)
    for(i in 1:length(Z)){
      F1[i]<-sum(Dec1$skews<=Z[i])/length(Dec1$skews)
      F2[i]<-sum(Dec2$skews<=Z[i])/length(Dec2$skews)
      F3[i]<-sum(Dec3$skews<=Z[i])/length(Dec3$skews)
    }
    Dec1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
    Dec2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
    Dec3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))

    ## Refit models
    gpd_boot_par <- GPD.fit(data_boot, gpd.quantile, optim.method)
    rate_boot_par <- rateparam.fit(data_boot, gpd.quantile, optim.method)
    extremalindex_boot_par <- extremalindex.fit(data_boot, EI.run.length, EI.quantile)

    period <- (max(data_boot$year)+1)-min(data_boot$year)
    base_yr <- min(data_boot$year)


    annmax.sealevel.dist <- function(z){
      F_k <- c()
      for(i in 1:period){
        F_m_s <- c()
        F_m_s[1] <- monmax.sealevel.dist(z,data,1,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[2] <- monmax.sealevel.dist(z,data,2,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[3] <- monmax.sealevel.dist(z,data,3,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[4] <- monmax.sealevel.dist(z,data,4,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[5] <- monmax.sealevel.dist(z,data,5,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[6] <- monmax.sealevel.dist(z,data,6,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[7] <- monmax.sealevel.dist(z,data,7,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[8] <- monmax.sealevel.dist(z,data,8,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[9] <- monmax.sealevel.dist(z,data,9,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[10] <- monmax.sealevel.dist(z,data,10,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[11] <- monmax.sealevel.dist(z,data,11,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)
        F_m_s[12] <- monmax.sealevel.dist(z,data,12,base_yr + (i-1), gpd_boot_par, rate_boot_par, extremalindex_boot_par, lin.int = EI_lin.interp, q=gpd.quantile, emp = empiricals)

        F_k[i] <- prod(F_m_s)
      }
      mean(F_k)-(1-p)
    }

    solver.interval.min <- quantile(data_boot$skews+data_boot$maxTide, 0.5, names=F)
    solver.interval.max <- max(data_boot$skews+data_boot$maxTide)+3
    while( sign(annmax.sealevel.dist(solver.interval.min)) == sign(annmax.sealevel.dist(solver.interval.max)) ){
      solver.interval.min <- solver.interval.min - 0.1
      solver.interval.max <- solver.interval.max + 0.1
    }
    RL[j] <- uniroot(annmax.sealevel.dist, c(solver.interval.min, solver.interval.max))$root
  }


  return(c(quantile(RL, (1-ci.prob)/2), quantile(RL, 0.5), quantile(RL, ci.prob+(1-ci.prob)/2)))
}

  
