#generalized distance function for n dimensions and unequal attribute weights
distfunct <- function(X,phi1,phi2,W,wgt) {
  
  n=dim(X)[1]
  D=matrix(NA,n,n)
  for (u in 1:n) {
    for (j in 1:n)  {
      distanz_allgemein = function(u,v,w,wgt){
        n = length(u)
        #Spannt den n-1-dimensionalen Indifferenzraum auf + Domimanzvektor
        B = matrix(0,n,n)
        #Calculate the indifference vectors
        for (i in 1:(n-1)){
          B[1,i] = -w[i+1]/w[1]
          B[i+1,i] = w[1]/w[1]
        }
        #calculate the doominancevector
        for (j in 1:n){
          B[j,n] = w[j]/w[1]
        }
        #Normierung der Basisvektoren auf die Länge 1:            
        for (k in 1:n){
          B[,k] = B[,k]/sqrt(sum((B[,k]^2)))
        }
        #Matrix zum berechnen des Abstandes = B
        #Den Abstand in den neuen Basisvektoren v1 bis v erhält man durch
        #Multiplikation mit P
        P = solve(B)
        #Berechnung des Abstandes
        abstand = u-v
        Distanz = P %*% abstand               
        #gewicht wgt einführen und abstand in der A-norm berechnen              
        A = diag(n)
        A[n,n] = wgt
        Distanz = t(Distanz) %*% A %*% Distanz               
        return(Distanz)            
      }
      D[u,j]= distanz_allgemein(X[u,],X[j,],W,wgt)
    } 
  } 
  D=.99*diag(n) - phi2 * exp(-phi1 * D^2)
  return(D)
} # end function: distfunct
