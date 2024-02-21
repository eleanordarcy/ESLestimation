get.empiricals <- function(data){
  Jan <- data %>% dplyr::filter(month==1)
  Feb <- data %>% dplyr::filter(month==2)
  Mar <- data %>% dplyr::filter(month==3)
  Apr <- data %>% dplyr::filter(month==4)
  May <- data %>% dplyr::filter(month==5)
  Jun <- data %>% dplyr::filter(month==6)
  Jul <- data %>% dplyr::filter(month==7)
  Aug <- data %>% dplyr::filter(month==8)
  Sep <- data %>% dplyr::filter(month==9)
  Oct <- data %>% dplyr::filter(month==10)
  Nov <- data %>% dplyr::filter(month==11)
  Dec <- data %>% dplyr::filter(month==12)
  
  Z <- seq(min(data$skews)-1,max(data$skews+1),by=0.001) #range of skews
  F1 <- F2 <- F3 <- c()
  ## Split skew surge empirical distribution into tidal bands associated with 0.33 and 0.67 quantiles
  q1 <- sort(Jan$maxTide)[0.333*length(Jan$maxTide)]
  q2 <- sort(Jan$maxTide)[0.667*length(Jan$maxTide)]
  Jan1 <- Jan %>% dplyr::filter( maxTide <= q1)
  Jan2 <- Jan %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Jan3 <- Jan %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Jan1$skews,3)<=Z[i])/length(Jan1$skews), 1/length(Jan1$skews))
    F2[i] <- max(sum(round(Jan2$skews,3)<=Z[i])/length(Jan2$skews), 1/length(Jan2$skews))
    F3[i] <- max(sum(round(Jan3$skews,3)<=Z[i])/length(Jan3$skews), 1/length(Jan3$skews))
  }
  Jan1_emp_skews <- lut(outputs=round(F1,3),inputs=round(Z,3))
  Jan2_emp_skews <- lut(outputs=round(F2,3),inputs=round(Z,3))
  Jan3_emp_skews <- lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1 <- F2<-F3<-c()
  q1 <- sort(Feb$maxTide)[0.333*length(Feb$maxTide)]
  q2 <- sort(Feb$maxTide)[0.667*length(Feb$maxTide)]
  Feb1 <- Feb %>% dplyr::filter( maxTide <= q1)
  Feb2 <- Feb %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Feb3 <- Feb %>% dplyr::filter( maxTide > q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Feb1$skews,3)<=Z[i])/length(Feb1$skews), 1/length(Feb1$skews))
    F2[i] <- max(sum(round(Feb2$skews,3)<=Z[i])/length(Feb2$skews), 1/length(Feb2$skews))
    F3[i] <- max(sum(round(Feb3$skews,3)<=Z[i])/length(Feb3$skews), 1/length(Feb3$skews))
  }
  Feb1_emp_skews <- lut(outputs=round(F1,3),inputs=round(Z,3))
  Feb2_emp_skews <- lut(outputs=round(F2,3),inputs=round(Z,3))
  Feb3_emp_skews <- lut(outputs=round(F3,3),inputs=round(Z,3))
  
  
  F1<-F2<-F3<-c()
  q1<-sort(Mar$maxTide)[0.333*length(Mar$maxTide)]
  q2<-sort(Mar$maxTide)[0.667*length(Mar$maxTide)]
  Mar1 <- Mar %>% dplyr::filter( maxTide <= q1)
  Mar2 <- Mar %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Mar3 <- Mar %>% dplyr::filter( maxTide > q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Mar1$skews,3)<=Z[i])/length(Mar1$skews), 1/length(Mar1$skews))
    F2[i] <- max(sum(round(Mar2$skews,3)<=Z[i])/length(Mar2$skews), 1/length(Mar2$skews))
    F3[i] <- max(sum(round(Mar3$skews,3)<=Z[i])/length(Mar3$skews), 1/length(Mar3$skews))
  }
  Mar1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Mar2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Mar3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Apr$maxTide)[0.333*length(Apr$maxTide)]
  q2<-sort(Apr$maxTide)[0.667*length(Apr$maxTide)]
  Apr1 <- Apr %>% dplyr::filter( maxTide <= q1)
  Apr2 <- Apr %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Apr3 <- Apr %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Apr1$skews,3)<=Z[i])/length(Apr1$skews), 1/length(Apr1$skews))
    F2[i] <- max(sum(round(Apr2$skews,3)<=Z[i])/length(Apr2$skews), 1/length(Apr2$skews))
    F3[i] <- max(sum(round(Apr3$skews,3)<=Z[i])/length(Apr3$skews), 1/length(Apr3$skews))
  }
  Apr1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Apr2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Apr3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(May$maxTide)[0.333*length(May$maxTide)]
  q2<-sort(May$maxTide)[0.667*length(May$maxTide)]
  May1 <- May %>% dplyr::filter( maxTide <= q1)
  May2 <- May %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  May3 <- May %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(May1$skews,3)<=Z[i])/length(May1$skews), 1/length(May1$skews))
    F2[i] <- max(sum(round(May2$skews,3)<=Z[i])/length(May2$skews), 1/length(May2$skews))
    F3[i] <- max(sum(round(May3$skews,3)<=Z[i])/length(May3$skews), 1/length(May3$skews))
  }
  May1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  May2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  May3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  
  F1<-F2<-F3<-c()
  q1<-sort(Jun$maxTide)[0.333*length(Jun$maxTide)]
  q2<-sort(Jun$maxTide)[0.667*length(Jun$maxTide)]
  Jun1 <- Jun %>% dplyr::filter( maxTide <= q1)
  Jun2 <- Jun %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Jun3 <- Jun %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Jun1$skews,3)<=Z[i])/length(Jun1$skews), 1/length(Jun1$skews))
    F2[i] <- max(sum(round(Jun2$skews,3)<=Z[i])/length(Jun2$skews), 1/length(Jun2$skews))
    F3[i] <- max(sum(round(Jun3$skews,3)<=Z[i])/length(Jun3$skews), 1/length(Jun3$skews))
  }
  Jun1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Jun2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Jun3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Jul$maxTide)[0.333*length(Jul$maxTide)]
  q2<-sort(Jul$maxTide)[0.667*length(Jul$maxTide)]
  Jul1 <- Jul %>% dplyr::filter( maxTide <= q1)
  Jul2 <- Jul %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Jul3 <- Jul %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Jul1$skews,3)<=Z[i])/length(Jul1$skews), 1/length(Jul1$skews))
    F2[i] <- max(sum(round(Jul2$skews,3)<=Z[i])/length(Jul2$skews), 1/length(Jul2$skews))
    F3[i] <- max(sum(round(Jul3$skews,3)<=Z[i])/length(Jul3$skews), 1/length(Jul3$skews))
  }
  Jul1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Jul2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Jul3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Aug$maxTide)[0.333*length(Aug$maxTide)]
  q2<-sort(Aug$maxTide)[0.667*length(Aug$maxTide)]
  Aug1 <- Aug %>% dplyr::filter( maxTide <= q1)
  Aug2 <- Aug %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Aug3 <- Aug %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Aug1$skews,3)<=Z[i])/length(Aug1$skews), 1/length(Aug1$skews))
    F2[i] <- max(sum(round(Aug2$skews,3)<=Z[i])/length(Aug2$skews), 1/length(Aug2$skews))
    F3[i] <- max(sum(round(Aug3$skews,3)<=Z[i])/length(Aug3$skews), 1/length(Aug3$skews))
  }
  Aug1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Aug2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Aug3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Sep$maxTide)[0.333*length(Sep$maxTide)]
  q2<-sort(Sep$maxTide)[0.667*length(Sep$maxTide)]
  Sep1 <- Sep %>% dplyr::filter( maxTide <= q1)
  Sep2 <- Sep %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Sep3 <- Sep %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Sep1$skews,3)<=Z[i])/length(Sep1$skews), 1/length(Sep1$skews))
    F2[i] <- max(sum(round(Sep2$skews,3)<=Z[i])/length(Sep2$skews), 1/length(Sep2$skews))
    F3[i] <- max(sum(round(Sep3$skews,3)<=Z[i])/length(Sep3$skews), 1/length(Sep3$skews))
  }
  Sep1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Sep2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Sep3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Oct$maxTide)[0.333*length(Oct$maxTide)]
  q2<-sort(Oct$maxTide)[0.667*length(Oct$maxTide)]
  Oct1 <- Oct %>% dplyr::filter( maxTide <= q1)
  Oct2 <- Oct %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Oct3 <- Oct %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Oct1$skews,3)<=Z[i])/length(Oct1$skews), 1/length(Oct1$skews))
    F2[i] <- max(sum(round(Oct2$skews,3)<=Z[i])/length(Oct2$skews), 1/length(Oct2$skews))
    F3[i] <- max(sum(round(Oct3$skews,3)<=Z[i])/length(Oct3$skews), 1/length(Oct3$skews))
  }
  Oct1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Oct2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Oct3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Nov$maxTide)[0.333*length(Nov$maxTide)]
  q2<-sort(Nov$maxTide)[0.667*length(Nov$maxTide)]
  Nov1 <- Nov %>% dplyr::filter( maxTide <= q1)
  Nov2 <- Nov %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Nov3 <- Nov %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Nov1$skews,3)<=Z[i])/length(Nov1$skews), 1/length(Nov1$skews))
    F2[i] <- max(sum(round(Nov2$skews,3)<=Z[i])/length(Nov2$skews), 1/length(Nov2$skews))
    F3[i] <- max(sum(round(Nov3$skews,3)<=Z[i])/length(Nov3$skews), 1/length(Nov3$skews))
  }
  Nov1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Nov2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Nov3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  F1<-F2<-F3<-c()
  q1<-sort(Dec$maxTide)[0.333*length(Dec$maxTide)]
  q2<-sort(Dec$maxTide)[0.667*length(Dec$maxTide)]
  Dec1 <- Dec %>% dplyr::filter( maxTide <= q1)
  Dec2 <- Dec %>% dplyr::filter( maxTide > q1 & maxTide <= q2)
  Dec3 <- Dec %>% dplyr::filter( maxTide >q2)
  for(i in 1:length(Z)){
    F1[i] <- max(sum(round(Dec1$skews,3)<=Z[i])/length(Dec1$skews), 1/length(Dec1$skews))
    F2[i] <- max(sum(round(Dec2$skews,3)<=Z[i])/length(Dec2$skews), 1/length(Dec2$skews))
    F3[i] <- max(sum(round(Dec3$skews,3)<=Z[i])/length(Dec3$skews), 1/length(Dec3$skews))
  }
  Dec1_emp_skews<-lut(outputs=round(F1,3),inputs=round(Z,3))
  Dec2_emp_skews<-lut(outputs=round(F2,3),inputs=round(Z,3))
  Dec3_emp_skews<-lut(outputs=round(F3,3),inputs=round(Z,3))
  
  ## Put into a named list
  empiricals <- list(Jan1_emp_skews, Jan2_emp_skews, Jan3_emp_skews,
                     Feb1_emp_skews, Feb2_emp_skews, Feb3_emp_skews,
                     Mar1_emp_skews, Mar2_emp_skews, Mar3_emp_skews,
                     Apr1_emp_skews, Apr2_emp_skews, Apr3_emp_skews,
                     May1_emp_skews, May2_emp_skews, May3_emp_skews,
                     Jun1_emp_skews, Jun2_emp_skews, Jun3_emp_skews,
                     Jul1_emp_skews, Jul2_emp_skews, Jul3_emp_skews,
                     Aug1_emp_skews, Aug2_emp_skews, Aug3_emp_skews,
                     Sep1_emp_skews, Sep2_emp_skews, Sep3_emp_skews,
                     Oct1_emp_skews, Oct2_emp_skews, Oct3_emp_skews,
                     Nov1_emp_skews, Nov2_emp_skews, Nov3_emp_skews,
                     Dec1_emp_skews, Dec2_emp_skews, Dec3_emp_skews)
  names(empiricals) <- c('Jan1_emp_skews', 'Jan2_emp_skews', 'Jan3_emp_skews',
                         'Feb1_emp_skews', 'Feb2_emp_skews', 'Feb3_emp_skews',
                         'Mar1_emp_skews', 'Mar2_emp_skews', 'Mar3_emp_skews',
                         'Apr1_emp_skews', 'Apr2_emp_skews', 'Apr3_emp_skews',
                         'May1_emp_skews', 'May2_emp_skews', 'May3_emp_skews',
                         'Jun1_emp_skews', 'Jun2_emp_skews', 'Jun3_emp_skews',
                         'Jul1_emp_skews', 'Jul2_emp_skews', 'Jul3_emp_skews',
                         'Aug1_emp_skews', 'Aug2_emp_skews', 'Aug3_emp_skews',
                         'Sep1_emp_skews', 'Sep2_emp_skews', 'Sep3_emp_skews',
                         'Oct1_emp_skews', 'Oct2_emp_skews', 'Oct3_emp_skews',
                         'Nov1_emp_skews', 'Nov2_emp_skews', 'Nov3_emp_skews',
                         'Dec1_emp_skews', 'Dec2_emp_skews', 'Dec3_emp_skews')
  return(empiricals)
}
