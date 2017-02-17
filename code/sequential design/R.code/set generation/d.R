###############################################
## Initial static stage                       #
## 1. Generate All profiles                   #
## 2. Select type of coding                   #
## 3. Select 5 sets based on prior betas      #
###############################################

#clear environment
rm(list=ls())

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
  
  if (length(betas)!=length(full_design)){return ("betas not correct")}
  
  #create vector with all Utilities
  N.row.design<-nrow(full_design)
  utilities<-as.vector(N.row.design)
  for (i in 1:N.row.design){utilities[i]<-sum(betas * full_design[i,])}  
  
  #vector with total utilities per set 
  total.U.set<-numeric(nrow(set.design))
  
  for (set in 1: nrow(set.design)){
    for (alt in 1: n_alts){
      prob.alt<-exp(utilities[set.design[set,alt]])
      total.U.set[set]<-total.U.set[set]+prob.alt
    }
  }
  
  #probabilities per alternative per set
  probs.per.set<-as.data.frame(matrix(nrow = nrow(set.design),ncol=n_alts))
  
  for (set in 1: nrow(set.design)){
    for (alt in 1: n_alts){
      probs.per.set[set,alt]<-exp(utilities[set.design[set,alt]])/total.U.set[set]
    }
  }
  
  return(probs.per.set)
  
}

# Information matrix for a choice set
Fish.inf<- function (set){
  
  prs.s<-as.numeric(probs.per.set[set,])
  #P1
  prs.diag<-diag(prs.s)
  
  #P1 - p1p1'
  M1<-(prs.s %*% t(prs.s))
  M2<-prs.diag - M1 
  
  #X'1(P1-p1p1')X1 = information matrix 
  X<-matrix(nrow=n_alts,ncol=length(betas))
  for(alt in 1:n_alts){ X[alt,]<-as.numeric(full_design[set.design[set,alt],])}
  
  M<-t(X) %*% M2  %*% X
  
  return(M)
  
}

# D-error for all choice sets
d_error<-function(){
  
  D<-as.vector(nrow(set.design))
  
  for (set in 1: nrow(set.design)){
    
    info<-Fish.inf(set)
    D[set]<-det(info+solve(COV))^(-1/length(betas))
  }
  
  return(D)
  
}


db_error<-function(){
  
  D<-as.vector(nrow(set.design))
  
  for (set in 1: nrow(set.design)){
    
    info<-Fish.inf(set)
    D[set]<-det(info+solve(COV))^(-1/length(betas))
  }
  return(D)
}






#Specify:-------------------------------------------------------
n_alts<-2                              # alternatives per set    
betas<-c(1.2,1.3)                      # Beta values             
beta.var<-c(0.8,1.1)                   # Beta variance 
COV<-diag(beta.var)
#---------------------------------------------------------------


#Functions
full_design<-full.profiles(levels_attributes=c(2,2),contrasts= "contr.sum")
all.choice.sets<-full.choice.sets()
probs.per.set<-probs()
d.error<-d_error()

#plot 
diff<-probs.per.set[,2]-probs.per.set[,1]
plot(diff,d.error)
text(diff, d.error, seq(1,length(d.error),1), pos=4, cex=0.7)





