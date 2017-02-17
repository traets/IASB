###
# Estimation expected human attention weights DFT (abad et al. 2014)
###

rm(list=ls())

#libraries
library(OpenMx)
library(mnormt)

######## SPECIFICATIONS #########
nAlts<-2        #Number alternatives
nAttr<-2        #Number attributes
Nd=500          #Number deliberation stepts (not specified in paper)

#responses (to reproduce table 2)
resp<-list(c(0,1),c(0,1),c(0,1),c(0,1),c(1,0),c(0,1),c(0,1),c(0,1),c(1,0))

#feedback matrix n x n 
s<-0.935             # not specified in paper (must be smaller than 1)
S<-diag(nAlts)
diag(S) <- s

#contrast matrix n x n
C = matrix(-1/(nAlts-1),nAlts,nAlts)
diag(C) = 1

#attribute matrix n x s (values of table 1)
att1A<-c(2.5, 3.74, 3.65, 3.12, 3.03, 3.26, 3.31, 3.01, 3.68)
att2A<-c(1.43, 1.4, 2.19, 1.4, 1.49, 1.43, 2.3, 2.11, 1.71)
att1B<-c(0.39, 0.9, 1.21, 1.46, 0.34, 1.71, 2.39, 2.33, 0.96)
att2B<-c(3.37, 3.65, 3.65, 3.17, 2.5, 3.15, 3.57, 4.04, 3.43)
MM<-cbind(att1A, att2A, att1B, att2B)

####### FUNCTIONS ##########
# L vector/matrix: contains -1 and 1, used to obtain mean difference in preference (see appendix B Roe et al. 2001) 
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

#mean preferences + Covariance matrix, returns probability
prob<-function(S,C,M,w,Nd){
  
  #set  
  I<-diag(nAlts)
  
  #mean
  Pref_mean <- solve(I-S) %*% (I-S^Nd)%*% C %*% M %*% w #(6) in Abad et al. I think it can be assumed that P(ti, 0)= 0 (some alternative specific preference?).
  #if t --> inf is assumed (I-S^Nd) also drops (as in Roe et al.) 
  
  #covariance valence (following appendix B of Roe et al.)
  covW<-  diag(w)-w%*%t(w)
  res<-diag(2); diag(res)<-0.1 #   Set as in Roe et al. Not sure how this is done in Abad et al.  
  
  cov_val<- C%*%M%*%covW%*%t(M)%*%t(C) + res 
  
  #covariance preference (Roe et al, appendix B)
  cov_pref<-matrix(data = 0, nrow=nAlts, ncol= nAlts)
  
  for(j in 0:(Nd-1)){
    A <- S^j%*%cov_val%*%t(S^j)
    cov_pref<- cov_pref + A
  }
  
  #probability alternative
  eta<-Pref_mean
  L<-L_fun(nAlts)
  
  probs <- apply(L,c(nAlts),proba, eta=eta, cov_pref=cov_pref)     # from berkowitsch 2014 code      
  
  return(probs)
  
}

###### ALGORITHM ########
#There are always 2 boundaries (A and B), the true one is in between. 
#Each boundary (slope) results in 2 attention weights, one for each attribute.
#The probabilities for choosing a particular alternative in a choice set (point in plane) can be calculated using the attention weights of each bound (thus twice).
#Using the EAW of the A boundary provides the lower probability for A choices and the upper boundary for B choices.
#Using the EAW of the B boundary provides the upper probability for A choices and the lower boundary for B choices.

#When a new choice is observed
# - falls between boundaries --> update one boundary dependent on which alternative is chosen.
# - falls in correct region (B beneath B bound, A above A bound) --> do nothing
# - conflict observation (B above A bound, A beneath B bound) --> compare likelihood of the conflicted observations 
# given we change the boundary with the likelihood of the present conflict observation given we do not change the boundary

#start slopes 
ba<--0
bb<--50 #minus infinity in abad et al. 
slopeA <- slopeB <- numeric(nrow(MM)) #vector for storing slopes (table 2 in abad) 

#plot the plane
plot(0, 0, ylim = c(-3,0), xlim = c(0,3))

#sequential updating 
for (i in 1:nrow(MM)){

  #choice set
  set<- MM[i,]
  M<-matrix(data = set, nrow = nAttr, ncol = nAlts, byrow = TRUE)
  
  #delta m1 & #delta m2
  delm1<-M[1,1]-M[2,1]
  delm2<-M[1,2]-M[2,2]
  
    # choice set between bounds --> update
    if (delm1*bb <= delm2 & delm2 <= delm1*ba){
      # A was chosen
      if (identical(resp[[i]], c(1,0))){
        ba<- delm2/delm1 
        A<-ba/(ba-1)
        wA<-c(A, 1-A)
    }
      # B was chosen
      if (identical(resp[[i]], c(0,1))){
        bb<- delm2/delm1
        B<-bb/(bb-1)
        wB<-c(B,1-B)
      }
    }
    
    # conflict: (A beneath B boundery, or B above A boundery)
    if ((delm2 < delm1*bb) && identical(resp[[i]], c(1,0)) || (delm1*ba < delm2) && identical(resp[[i]], c(0,1))){
      
      # A was chosen
      if (identical(resp[[i]], c(1,0))){
      P1<-prob(S=S, C=C, M=M, w=wB, Nd=Nd)[1]
        
        # quadrant 2
        # working on 
        
        # quadrant 4
        banew<- delm2/delm1
        newA<-banew/(banew-1)
        wAnew<-c(newA,1-newA)
        
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
            if(DELm2 > DELm1* banew){
              
              #calculate probability
              pi<-prob(S=S, C=C, M=M, w=wAnew, Nd=Nd)[2]
            }else{pi<-1}
            
            #multiply all probs (likelihood of conflict observations)
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
        # working on 
      
        # quadrant 4
        bbnew<- delm2/delm1
        newB<-bbnew/(bbnew-1)
        wBnew<-c(newB,1-newB)
        
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
              if(DELm2< DELm1* bbnew){
              pi<-prob(S=S, C=C, M=M, w=wBnew, Nd=Nd)[2]
              }else{pi<-1}
              
            #multiply all probs 
            P2<-P2*pi
          }
        }
        
        #compare
        if (P2 > P1){ bb<-bbnew } #change occured
        
    }
        
    }
  
  #save estimated slopes (table 2) 
  slopeA[i]<-ba
  slopeB[i]<-bb
  
  #plot
  #label in plot
  if (identical(resp[[i]], c(1,0))){
    label = LETTERS[1]
  }else{label = LETTERS[2]}
  
  points(delm1, delm2, pch=label)
  abline(a=0, b=ba,col="blue")
  abline(a=0,b=bb,col="red")

}

### Table 2
table2<-cbind(slopeA, slopeB)
table2

### Table 3 probability for choosing A (using estimates after t8)
B_t8<-c(0.374359, 0.625641)
A_t8<-c(0.272973, 0.727027)

#matrix to store probabilities
probis<-matrix(data = NA, nrow= 8, ncol=2)

for (i in 1:8){
  
  #choice set
  set<- MM[i,]
  M<-matrix(data = set, nrow = nAttr, ncol = nAlts, byrow = TRUE)
  
  #probability choice = A
  probis[i,1]<- probA<-prob(S=S, C=C, M=M, w=A_t8, Nd=Nd)[1]
  probis[i,2]<- probB<-prob(S=S, C=C, M=M, w=B_t8, Nd=Nd)[1]

}

options(scipen=999)
probis  
# Er zitten nog kleine verschillen op die waarschijnlijk te maken hebben met het feit dat Nd en de S diagonaal 
# niet hetzelfde zijn (waren niet gegeven).   



