
# gegeben Matrix S und Matrix Phi (beide 3x3)
#
# Aufstellen der Matrix A, S * Omega * S' = A * Omega
#
#

omegaNoLoop = function(S,Phi,nAlt) {
    
    #Step 1 : Aufstellen der Matrix I-A (9x9)
    A = matrix(NA,nAlt^2,nAlt^2)
      
    for (i in 1:nAlt){
        for (j in 1:nAlt) {
            A[(1+(i-1)*nAlt):(nAlt*i),(1+(j-1)*nAlt):(nAlt*j)] = S[i,j] * S
        }
    }
    
    I = diag(nAlt^2)
    A = (I-A)
    
    #Schritt 3 : Aufstellen vom Vektor Psi (9x1)
    
    Phi = t(Phi)
    dim(Phi) = c(nAlt^2,1)
        
    #Schritt 4 : Lösen des Gleichungssystems, Transformation auf 3x3 Matrix
    
#     if(abs(det(A))<10E-15) {
    if(abs(det(A))<10E-17) {
        #print("A is singulary")
        OmegaDirect = matrix(NA,nAlt,nAlt)
    } else {
        Omega = try(solve(A,Phi), 
                    silent = T)
        OmegaDirect = t(Omega)
        dim(OmegaDirect) = c(nAlt,nAlt)    
    }

    return(OmegaDirect)
       
}
