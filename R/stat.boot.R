stat.boot <- function(X,mean.block.size){
  N=length(X)
  
  #Generate random block sizes
  b=rgeom(N,1/(mean.block.size))+1
  
  sum=cumsum(b)
  
  b=b[1:min(which(sum>=N))]
  
  #Find starting indices
  inds=sample(1:N,length(b))
  all_inds=c()
  for(i in 1:length(b)){
    
    block_inds=inds[i]:(inds[i]+b[i]-1)
    
    #Wrap around indices
    if(sum(block_inds>N)>0){
      block_inds[block_inds>N]=1:sum(block_inds>N)
    }
    all_inds=c(all_inds,block_inds)
    
  }
  
  #Cut down to N only
  all_inds=all_inds[1:N]
  
  boot=X[all_inds]
  return(boot)
}
