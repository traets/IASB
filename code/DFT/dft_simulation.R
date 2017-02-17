###
# simulation
###
rm(list=ls())

#libraries
library(OpenMx)
library(mnormt)
library(mlogit)

#set.seed(0105)

######## SPECIFICATIONS #########
nAlts<-2        #Number alternatives
nAttr<-2        #Number attributes
Nd=500          #Number deliberation stepts (not specified in paper)
N<-50           #Number choice sets

#true w
w_true<-c(0.36, 0.64)

#feedback matrix n x n 
s<-0.935             # not specified in paper (must be smaller than 1)
S<-diag(nAlts)
diag(S) <- s

#contrast matrix n x n
C = matrix(-1/(nAlts-1),nAlts,nAlts)
diag(C) = 1

#respons vector
resp<-list()

##design
att1A<-runif(N, 0, 4)
att2A<-runif(N, 0, 4)
att1B<-runif(N, 0, 4)
att2B<-runif(N, 0, 4)
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

#respond function
respond<-function(S,C,M,w,Nd){
  
  #delta m1 & #delta m2
  delm1<-M[1,1]-M[2,1]
  delm2<-M[1,2]-M[2,2]
  
  ## based on quadrant transform or decide
  #quadrant 1 --> A
  if (delm1 > 0 && delm2 > 0){resp<-c(1,0)} 
  #quadrant 2 --> transform to 4
  else if (delm1 < 0 && delm2 > 0){
    delm2_new=-abs(delm1)
    delm1_new=delm2
    
    M_new<-matrix(data = NA, nrow = nAttr, ncol = nAlts, byrow = TRUE)
    M_new[,1]<-c(1+delm1_new, 1)
    M_new[,2]<-c(1+delm2_new, 1)
    
    pA<-prob(S=S, C=C, M=M_new, w=w, Nd=Nd)[1]
   
    A<-runif(1,0,1)
    if ( A < pA ){resp <- c(1,0) }else{resp<-c(0,1)}
    
  }
  #quadrant 3 --> B
  else if (delm1 < 0 && delm2 < 0){resp<-c(0,1)}
  #quadrant 4
  if (delm1 > 0 && delm2 < 0){
  
  pA<-prob(S=S, C=C, M=M, w=w, Nd=Nd)[1]
  A<-runif(1,0,1)
  if (A < pA ){resp <- c(1,0) }else{resp<-c(0,1)}
  }
  
  return(resp)
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
plot(0, 0, ylim = c(-3,3), xlim = c(-3,3))
abline(h=0)
abline(v=0)
abline(a=0, b=-0.56,col="green") #true DBL

#sequential updating 
for (i in 1:nrow(MM)){
  
  #Sys.sleep(1.5)
  
  #choice set
  set<- MM[i,]
  M<-matrix(data = set, nrow = nAttr, ncol = nAlts, byrow = TRUE)
  
  #response with true weights
  resp[[i]]<-respond(S=S, C=C, M=M, w = w_true, Nd=Nd)
  
  ###Estimate
  
  #delta m1 & #delta m2
  delm1<-M[1,1]-M[2,1]
  delm2<-M[1,2]-M[2,2]
  
  
  #Quadrant 2
  if (delm1 < 0 && delm2 > 0){
    
    x=-abs(delm1)
    y=delm2
    
    delm1<-y
    delm2<-x
    points(delm1, delm2, col="green")
  }
  
  #Quadrant 4
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
      P1<-prob(S=S, C=C, M=M, w=wA, Nd=Nd)[2]
      
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


#estimated boundary between
wA
wB
w_true

# #Design matrix + responsevector
# Y<-numeric(N)
# for (i in 1:N){
#   if(identical(resp[[i]], c(1,0)))
#   {Y[i]<-1}else {Y[i]<-2}
# } 
# 
# data<-as.data.frame(cbind(MM,Y))
# colnames(data)<- c("attr1.A", "attr2.A", "attr1.B", "attr2.B", "Y")
# 
# #Estimate multinomial logit 
# data_mnl<-mlogit.data(data, choice="Y", shape="wide", varying =1:4, alt.levels = 1:2, sep = "."  )
# 
# MNL_1 <- mlogit(Y ~ attr1 + attr2|-1, data_mnl, method="bfgs", print.level=0)
# coefs<-as.numeric(MNL_1$coefficients)
# 
# weights<-coefs/sum(coefs)
# weights
# wA
# wB
# w_true




