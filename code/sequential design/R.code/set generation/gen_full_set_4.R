###############################################
## Initial static stage
## 1. Generate All profiles
## 2. Select type of coding (work in progr)
## 3. Select 5 best sets based on prior betas
###############################################

#clear environment
rm(list=ls())

#1.state the nr of attributes and their levels (with a vector of all the nr of levels.) 
##Function generating ALL profiles.----------------------
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

full_design<-full.profiles(levels_attributes=c(3,4,4,5),contrasts= "contr.sum")



### 2. Specify prior. Select N most efficient sets, based on D-efficiency  

#specify 
n_alts<-3                       #n_alts in one choice set
betas<-c(0.805,-0.080,-0.603)   #prior = vector of part-worths 
fun<-function (x){return(1:x)}
e<-exp(1) 
N.row.design<-nrow(full_design)


#Create all sets 
sets.1<-rep(N.row.design,n_alts)
sets.2<-as.list(sets.1)
sets.3<-lapply(X = sets.2, fun)
set.design<-as.data.frame(expand.grid(sets.3))

#generate probalilities -----------------------------------------
#create vector with all Utilities
utilities<-as.vector(N.row.design)
for (i in 1:N.row.design){utilities[i]<-sum(betas*full_design[i,])}  

#vector with total utilities per set 
total.U.set<-numeric(nrow(set.design))

for (set in 1: nrow(set.design)){
  for (alt in 1: n_alts){
    prob.alt<-e^utilities[set.design[set,alt]]
    total.U.set[set]<-total.U.set[set]+prob.alt
  }
}

#probabilities per alternative per set
probs.per.set<-as.data.frame(matrix(nrow = nrow(set.design),ncol=n_alts))

for (set in 1: nrow(set.design)){
  for (alt in 1: n_alts){
    probs.per.set[set,alt]<-e^utilities[set.design[set,alt]]/total.U.set[set]
  }
}

#D-error for a set---------------------------
#p1
setteller<-7
prs.s<-as.numeric(probs.per.set[setteller,])
#P1
prs.diag<-diag(prs.s)

#P1 - p1p1'
mat1<-(prs.s %*% t(prs.s))
mat2<-prs.diag - mat1 

#X'1(P1-p1p1')X1
X<-matrix(nrow=n_alts,ncol=length(betas))
for(alt in 1:n_alts){ X[alt,]<-as.numeric(full_design[set.design[setteller,alt],])}

mat3<-t(X) %*% mat2 
mat4<-mat3%*%X
mat4
#-----------------------------------------



(18*18)-18



