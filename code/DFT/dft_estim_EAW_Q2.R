

###
# Estimation expected human attention weights DFT
###
rm(list=ls())

#library
library(LaplacesDemon)
library(OpenMx)
library(mnormt)

######## SPECIFICATIONS #########
nAlts<-2
nAttr<-2
Nd=500

#responses
resp<-list(c(0,1),c(0,1),c(0,1),c(0,1),c(1,0),c(0,1),c(0,1),c(0,1),c(1,0))

#feedback matrix n x n 
s<-0.88
S<-diag(nAlts)
diag(S) <- s
#contrast matrix n x n
C = matrix(-1/(nAlts-1),nAlts,nAlts)
diag(C) = 1

#attribute matrix n x s
att1A<-c(2.5, 3.74, 3.65, 3.12, 3.03, 3.26, 3.31, 3.01, 3.68)
att2A<-c(1.43, 1.4, 2.19, 1.4, 1.49, 1.43, 2.3, 2.11, 1.71)
att1B<-c(0.39, 0.9, 1.21, 1.46, 0.34, 1.71, 2.39, 2.33, 0.96)
att2B<-c(3.37, 3.65, 3.65, 3.17, 2.5, 3.15, 3.57, 4.04, 3.43)
MM<-cbind(att1A, att1B, att2A, att2B)

####### FUNCTIONS ##########
# L vector/matrix
L_fun<-function(nAlts){
  
  L=array(NA,c(nAlts-1,nAlts,nAlts))
  counter=1   
  for(k in 1:nAlts){
    for(i in 1:(nAlts-1)){
      for(j in 1:nAlts){
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
  return(L)
}

# choice probability
proba<-function(L, eta, cov_pref){
  Leta = L%*%eta # mean difference in preference
  Lomega = L%*%cov_pref%*%t(L) # varance-covariance matrix of the mean difference in preference
  p = pmnorm(x=Leta,varcov=Lomega)
  
  return(p)
}
prob<-function(S,C,M,w,Nd){
  
  #set  
  I<-diag(nAlts)
  
  #mean
  Pref_mean <- solve(I-S) %*% (I-S^Nd)%*% C %*% M %*% w 
  
  #covariance valence(t)
  covW<-  diag(w)-w%*%t(w)
  cov_val<- C%*%M%*%covW%*%t(M)%*%t(C)
  
  #cov preference
  cov_pref<-matrix(data = 0, nrow=nAlts, ncol= nAlts)
  
  for(j in 0:(Nd-1)){
    A <- S^j%*%cov_val%*%t(S^j)
    cov_pref<- cov_pref + A
  }
  
  #probability alternative
  eta<-Pref_mean
  L<-L_fun(nAlts)
  
  probs <- apply(L,c(nAlts),proba, eta=eta, cov_pref=cov_pref)          
  
  return(probs)
  
}

###### ALGORITHM ########

#start slopes
ba<--0
bb<--50

#plot
plot(0, 0, ylim = c(0,3), xlim = c(-3,0))

for (i in 1:nrow(MM)){
i=1
  
  #choice set
  set<- MM[i,]
  M<-matrix(data = set, nrow = nAttr, ncol = nAlts, byrow = TRUE)
  M
  
  #delta m1 & #delta m2
  delm1<-M[1,1]-M[2,1]
  delm2<-M[1,2]-M[2,2]
  
  # choice set between bounds --> update
  if (delm2*ba <= delm1 & delm1 <= delm2*bb){
    # A was chosen
    if (identical(resp[[i]], c(1,0))){
      ba<- delm2/delm1 
      p1A<-ba/(ba-1)
      p_A<-c(p1A,1-p1A)
    }
    # B was chosen
    if (identical(resp[[i]], c(0,1))){
      bb<- delm2/delm1
      p1B<-bb/(bb-1)
      p_B<-c(p1B,1-p1B)
    }
  }
  
  # conflict: (A beneath B boundery, or B above A boundery)
  if ((delm2 < delm1*bb) && identical(resp[[i]], c(1,0)) || (delm2*ba < delm1) && identical(resp[[i]], c(0,1))){
    
    # A was chosen
    if (identical(resp[[i]], c(1,0))){
      P1<-prob(S=S, C=C, M=M, w=p_B, Nd=Nd)[1]
      
      # quadrant 2
      
      # QUADRANT 4
      banew<- delm2/delm1
      pnewA<-banew/(banew-1)
      p_Anew<-c(pnewA,1-pnewA)
      
      #probability inconsistent set of B responses according to new ba 
      P2<-1
      for (k in 1:i){
        
        
        #for all previous sets with B resonses
        if(identical(resp[[k]], c(0,1))){
          
          setk<- MM[k,]
          M<-matrix(data = setk, nrow = nAttr, ncol = nAlts, byrow = TRUE)
          
          DELm1<-M[1,1]-M[2,1]
          DELm2<-M[1,2]-M[2,2]
          
          #that fall above the new ba slope
          if(DELm1 > DELm2* banew){
            pi<-prob(S=S, C=C, M=M, w=p_Anew, Nd=Nd)[2]
          }else{pi<-1}
          
          #multiply all probs 
          P2<-P2*pi
        }
      }
      
      #compare
      if (P2 > P1){ ba<-banew } #change occured
      
    }
    
    # B was chosen
    if (identical(resp[[i]], c(0,1))){
      P1<-prob(S=S, C=C, M=M, w=p_A, Nd=Nd)[2]
      
      # quadrant 2
      
      # quadrant 4
      bbnew<- delm2/delm1
      pnewB<-bbnew/(bbnew-1)
      p_Bnew<-c(pnewB,1-pnewB)
      
      #probability inconsistent set of A responses according to new bb 
      P2<-1
      for (k in 1:i){
        
        #for all previous sets with A resonses
        if(identical(resp[[k]], c(1,0))){
          
          setk<- MM[k,]
          M<-matrix(data = setk, nrow = nAttr, ncol = nAlts, byrow = TRUE)
          
          DELm1<-M[1,1]-M[2,1]
          DELm2<-M[1,2]-M[2,2]
          
          #that fall beneath the new BB slope
          if(DELm1< DELm2* bbnew){
            pi<-prob(S=S, C=C, M=M, w=p_Bnew, Nd=Nd)[2]
          }else{pi<-1}
          
          #multiply all probs 
          P2<-P2*pi
        }
      }
      
      #compare
      if (P2 > P1){ bb<-bbnew } #change occured
      
    }
    
  }
  
  
  ###save probs
  p1A<-ba/(ba-1)
  p_A<-c(p1A,1-p1A)
  probA<-prob(S=S, C=C, M=M, w=p_A, Nd=Nd)
  
  
  p1B<-bb/(bb-1)
  p_B<-c(p1B,1-p1B)
  probB<-prob(S=S, C=C, M=M, w=p_B, Nd=Nd)
  
  
  ###plot
  if (identical(resp[[i]], c(1,0))){
    label = LETTERS[1]
  }else{label = LETTERS[2]}
  
  points(delm2, delm1, pch=label)
  abline(a=0, b=ba,col="blue")
  abline(a=0,b=bb,col="red")
  
}


### Table 3 ###
#matrix to store probabilities
# probis<-matrix(data = NA, nrow= nrow(MM), ncol=2)
# 
# #probability for choosing A (using estimates after t8)
# B_t8<-c(0.374359, 0.625641)
# A_t8<-c(0.272973, 0.727027)
# 
# for (i in 1:8){
#   
#   #choice set
#   set<- MM[i,]
#   M<-matrix(data = set, nrow = nAttr, ncol = nAlts, byrow = TRUE)
#   
#   #probability choice = A
#   probis[i,1]<- probA<-prob(S=S, C=C, M=M, w=A_t8, Nd=Nd)[1]
#   probis[i,2]<- probB<-prob(S=S, C=C, M=M, w=B_t8, Nd=Nd)[1]
# }
# 
# probis
# 


