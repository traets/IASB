#############################
# calculate Db- error with  #
# swapping algorithm.       #
#############################

#clear environment
rm(list=ls())

#libraries
library(dplyr)

### 1. State the nr of attributes and their levels (with a vector of all the nr of levels.) 

#Function generating ALL profiles.----------------------
full.profiles<-function (levels_attributes,contrasts) {
  
  #as list
  levels_attributes<-as.list(levels_attributes)
  
  #Make list with all atribute levels
  fun<-function (x){return(1:x)}
  levels_attributes<-lapply(X = levels_attributes, fun)
  
  #generate every combination
  full_design<-as.data.frame(expand.grid(levels_attributes))
  
  #create factors:
  col_names <- names(full_design)
  full_design[,col_names] <- lapply(full_design[,col_names] , factor)
  
  #contrast list
  con <- list()
  for(i in 1:length(levels_attributes)){
    name <- paste('Var',i,sep='')
    con[name] <- contrasts
  }
  
  #coding
  full_D<-as.data.frame(model.matrix(~ ., full_design, contrasts = con))
  #delete intercept
  full_D<-full_D[,-1]
  
  return(full_D)
  
}

# Create all possible choice sets 
full.choice.sets<- function(){
  
  fun<-function (x){return(1:x)}
  N.row.design<-nrow(full_design)
  
  sets.1<-rep(N.row.design,n_alts)
  sets.2<-as.list(sets.1)
  sets.3<-lapply(X = sets.2, fun)
  set.design<-as.data.frame(expand.grid(sets.3))
  
  #to global env
  assign("set.design", set.design, envir=globalenv())
  
  #return
  return (set.design)
}

# Generate probalilities 
probs<-function (){

  #create vector with all Utilities
  utilities<-as.vector(nrow(S1))
  for (i in 1:nrow(S1)){utilities[i]<-exp(sum(betas * S1[i,]))}  
  
  #vector with total utilities per set 
  total.U.set<-vector(mode="numeric", length=nrow(S1))
  setnr<-seq(1,nrow(S1),n_alts)
  
  for (set in setnr){
   
    for(alts in 1:n_alts){
  
      U.alt<-utilities[set+(alts-1)]
      total.U.set[set]<-(total.U.set[set]+U.alt)
      total.U.set[set+(alts-1)]<-total.U.set[set]
   }
  }
  
  #probabilities per alternative per set
  probs.per.set<-utilities/total.U.set
  
  return(probs.per.set)
}

# Information matrix for a choice set
Fish.inf<- function (){
  
  setnr<-seq(1,nrow(S1),n_alts)
  prs<-probs()
  prs.s<-split(prs,setnr)
  prs.s  
  set=2
  for (set in setnr){
  #P1
  prs.diag<-diag(prs.s[set])
  
  #P1 - p1p1'
  M1<-(prs.s %*% t(prs.s))
  M2<-prs.diag - M1 
  
  #X'1(P1-p1p1')X1 = information matrix 
  X<-matrix(nrow=n_alts,ncol=length(betas))
  for(alt in 1:n_alts){ X[alt,]<-as.numeric(full_design[set.design[set,alt],])}
  
  M[setnr]<-t(X) %*% M2  %*% X
}
  
}

# D-error for all choice sets
d.error<-function(){
  
  D<-as.vector(nrow(n_init.sets))
  
  for (set in 1: nrow()){
    
    info<-Fish.inf(set)
    D[set]<-det(info+solve(COV))^(-1/length(betas))
  }
  
  return(D)
  
}

#Specify:
n_alts<-2                                 # alternatives per set
betas<-c(-0.5,1,0.5,0.8,1.2)                      # Beta values 
beta.var<-c(2,2,2,1,1)                        # Beta variance 
COV<-diag(beta.var)
n_init.sets<-5                           #nr of initial choice sets per pp
N_resp<-10                               #nr of respondents 


#1. generate all possible profiles
full_profiles<-full.profiles(levels_attributes=c(3,3,2),contrasts= "contr.sum")

#2. create a random subdesign containing n_init.sets from full_design.
S1<-matrix(nrow=(n_alts*n_init.sets),ncol=length(betas))
S1<-sample_n(full_profiles, 10)

#3. optimize d-error of subdesign by swapping
    #3.1. calculate d-error for subdesign S1
      #probability vector for Subdesign 
      probs.sub<-probs()



#4. create next subdesign conditionally on previous one
#5. repeat untill N_resp subdesigns 




