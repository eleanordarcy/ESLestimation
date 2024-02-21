emp.dist <- function(data, value, obs.tide, mnth, emp ){
  q1 <- quantile((data %>% dplyr::filter(month == mnth))$maxTide, 0.333, names=F) 
  q2 <- quantile((data %>% dplyr::filter(month == mnth))$maxTide, 0.667, names=F) 
  
  if(mnth==1){
    if(obs.tide<=q1){
      return(emp$Jan1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Jan2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Jan3_emp_skews(round(value,3)))
    }
  } else if(mnth == 2){
    if(obs.tide<=q1){
      return(emp$Feb1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Feb2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Feb3_emp_skews(round(value,3)))
    }
  } else if(mnth == 3){
    if(obs.tide<=q1){
      return(emp$Mar1_emp_skews(round(value,3)))
    } 
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Mar2_emp_skews(round(value,3)))
    } 
    else{
      return(emp$Mar3_emp_skews(round(value,3)))
    }
  } else if(mnth == 4){
    if(obs.tide<=q1){
      return(emp$Apr1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Apr2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Apr3_emp_skews(round(value,3)))
    }
  } else if(mnth == 5){
    if(obs.tide<=q1){
      return(emp$May1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$May2_emp_skews(round(value,3)))
    }
    else{
      return(emp$May3_emp_skews(round(value,3)))
    }
  } else if(mnth == 6){
    if(obs.tide<=q1){
      return(emp$Jun1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Jun2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Jun3_emp_skews(round(value,3)))
    }
  } else if(mnth == 7){
    if(obs.tide<=q1){
      return(emp$Jul1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Jul2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Jul3_emp_skews(round(value,3)))
    }
  } else if(mnth == 8){
    if(obs.tide<=q1){
      return(emp$Aug1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Aug2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Aug3_emp_skews(round(value,3)))
    }
  } else if(mnth == 9){
    if(obs.tide<=q1){
      return(emp$Sep1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Sep2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Sep3_emp_skews(round(value,3)))
    }
  } else if(mnth == 10){
    if(obs.tide<=q1){
      return(emp$Oct1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Oct2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Oct3_emp_skews(round(value,3)))
    }
  } else if(mnth == 11){
    if(obs.tide<=q1){
      return(emp$Nov1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Nov2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Nov3_emp_skews(round(value,3)))
    }
  } else if(mnth == 12){
    if(obs.tide<=q1){
      return(emp$Dec1_emp_skews(round(value,3)))
    }
    else if(obs.tide > q1 & obs.tide <= q2){
      return(emp$Dec2_emp_skews(round(value,3)))
    }
    else{
      return(emp$Dec3_emp_skews(round(value,3)))
    }
  }
}
