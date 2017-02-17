####################################
## Initial static stage
## 1. Generate All profiles
## 2. Select type of coding (work in progr)
## 3. Select 5 best sets based on prior 
####################################

#clear environment
rm(list=ls())

#1.state the nr of attributes and their levels (with a vector of all the nr of levels.) 
##Function generating ALL profiles.----------------------
full.profiles<-function (levels_attributes) {

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
 
  
  #coded
  full_D<-as.data.frame(model.matrix(~ ., full_design, contrasts = list(Var1="contr.sum", Var2="contr.sum")))
  #delete intercept
  full_D<-full_D[,-1]

  
  return(full_D)

}


full_design<-full.profiles(levels_attributes=c(3,2))



#2. Specify prior. Select N most efficient sets, based on D-efficiency  

#specify 
n_alts<-2                       #n_alts in one choice set
betas<-c(0.805,-0.080,-0.603)   #prior = vector of part-worths 

e<-exp(1) 
N.row.design<-nrow(full_design)

utilities<-as.vector(N.row.design)
#prob_set<-data.frame(matrix(nrow=N.row.design^n_alts,ncol= n_alts))
  
# Create vector with all Utilities  
for (i in 1:N.row.design){utilities[i]<-sum(betas*full_design[i,])}  

#probabilities for every design 
sets.1<-rep(N.row.design,n_alts)
sets.2<-as.list(sets.2)

#Make list with all atribute levels
fun<-function (x){return(1:x)}
sets.3<-lapply(X = sets.2, fun)

#generate every combination
set.design<-as.data.frame(expand.grid(sets.3))

#generate probalilities 


