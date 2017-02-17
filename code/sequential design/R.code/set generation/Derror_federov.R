#######################################################################
## Calculate D-efficiency of a choice set given some specifications   #
## 1. Specify betas, nr of alts/choice set, nr attributes, levels     #
## 2. calculate D-error one set at a time                             #
## 3. change alternatives of set to improve D.error                   #
## 4. Store nr of desired initial sets                                #
#######################################################################


#clear environment
rm(list=ls())

### 1. State the nr of attributes and their levels (with a vector of all the nr of levels.) 

#Function generating ALL profiles.
full.profiles<-function (contrasts) {
  
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
#select one start choice set
random.choice.set<-function (){
  
  choice.set<-matrix(nrow = n_alts,ncol = length(all_profiles))
  
  for (i in 1:n_alts){
    r<-round(runif(1,1,nrow(all_profiles)))
    choice.set[i,]<-as.numeric(all_profiles[r,]) }
  
  return(choice.set)
}
# Generate probalilities 
probs<-function (){
  
  if (length(betas)!=ncol(choice_set)){stop ("nr of betas not correct")}
  
  #create vector with all Utilities
  utilities<-as.vector(n_alts)
  for (i in 1:n_alts){utilities[i]<-sum(betas * choice_set[i,])}  
  
  #vector with total utilities per set 
  total.U.set<-sum(exp(utilities))
  
  #probs for alts 
  probs<-as.numeric(n_alts)
  for (a in 1:n_alts){probs[a]<-exp(utilities[a])/total.U.set}
  
  return(probs)
  
}
# d.error for a choice set
d.error<- function (){
  
  probs.set<-probs()
  prs.s<-probs.set
  #P1
  prs.diag<-diag(prs.s)
  
  #P1 - p1p1'
  M1<-(prs.s %*% t(prs.s))
  M2<-prs.diag - M1 
  
  #X'1(P1-p1p1')X1 = information matrix 
  X<-choice_set
  M<-t(X) %*% M2  %*% X
  
  #d-error
  D<-det(M+solve(COV))^(-1/length(betas))
  
}


#######################################################################################################

#Specify:-------------------------------------------------------
levels_attributes=c(3,3,3)
n_alts<-2                                          # alternatives per set    
betas<-c(0.2, 1.2, 0.3, 1.8, 0.8, 0.5)              # Beta values             
beta.var<-c(2, 1, 1.5, 0.8, 1.2, 1.9)                   # Beta variance 
COV<-diag(beta.var)
n_initial.sets<-5                                  #nr of best sets to store
best.error<-4
best.error2<-2
best.error3<-1.5
best.error4<-1.2
best.error5<-1

K<-1000
#---------------------------------------------------------------

#1.calculate d-error given betas/ choice set/ COV
#generate full_design
all_profiles<-full.profiles(contrasts= "contr.sum")
#Store N best choice sets
sets<-data.frame(matrix(nrow = n_alts*n_initial.sets, ncol=length(betas)))


for (i in 1:K)
{
#draw random initial set
choice_set<-random.choice.set()

#d.error
D_error<-d.error()

if (D_error< best.error){
 for (i in 1:n_alts){
 sets[i,]<-as.numeric((choice_set[i,]))
 }
  best.error<-D_error
  
}

else if (D_error< best.error2){
  for (i in 1:n_alts){
    sets[i+n_alts,]<-as.numeric((choice_set[i,]))
  }
  best.error2<-D_error
}

else if (D_error< best.error3){
  for (i in 1:n_alts){
    sets[i+(2*n_alts),]<-as.numeric((choice_set[i,]))
  }
  best.error3<-D_error
}

else if (D_error< best.error4){
  for (i in 1:n_alts){
    sets[i+(3*n_alts),]<-as.numeric((choice_set[i,]))
  }
  best.error4<-D_error
}

else if (D_error< best.error5){
  for (i in 1:n_alts){
    sets[i+(4*n_alts),]<-as.numeric((choice_set[i,]))
  }
  best.error5<-D_error
}

print(choice_set)
}


sets
best.error
