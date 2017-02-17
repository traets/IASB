#create the C-Matrix and the L-array, based on the number of alternatives


C = matrix(-1/(nAlt-1),nAlt,nAlt)
diag(C) = 1

counter=1    
L=array(NA,c(nAlt-1,nAlt,nAlt))
for(k in 1:nAlt){
  for(i in 1:(nAlt-1)){
    for(j in 1:nAlt){
      if(counter==k && counter==j){
        counter=counter+1
      }
      if(j==k){
        L[i,j,k]=1
      } else {
        if (j==counter){
          L[i,j,k]=-1
        } else {
          L[i,j,k]=0
        }
      }  
    }
    counter=counter+1
  }
  counter=1
}
