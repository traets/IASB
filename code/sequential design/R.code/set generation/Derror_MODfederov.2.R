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

#generate start design 
start.design<-function (){
  
  start<-data.frame(matrix(nrow=n_sets*n_alts, ncol=length(betas)))
  
  for(i in 1:nrow(start)){
    R<-round(runif(1, 1, nrow(candidate_profiles)))
    start[i,]<-as.numeric(candidate_profiles[R,])
  }
  
  
  set<-seq(1, n_sets, 1)
  start$set<-rep(set, each=n_alts)
  
  alt<-seq(1, n_alts, 1)
  start$alt<-rep(alt, n_sets)
  
  #return 
  return(start)
}

# d.error for a choice set
sum.information<- function () {
  
  for ( c_set in 1:n_sets){
    
    
    drops <- c("set","alt")
    choice_set<-start[start$set==c_set,!(names(start) %in% drops)]
  
    ####information
    #create vector with all Utilities
    utilities<-as.vector(n_alts)
    for (i in 1:n_alts){utilities[i]<-sum(betas * choice_set[i,])}  
    
    #vector with total utilities per set 
    total.U.set<-sum(exp(utilities))
    
    #probs for alts 
    prs.s<-as.numeric(n_alts)
    for (a in 1:n_alts){prs.s[a]<-exp(utilities[a])/total.U.set}
  
    #P1
    prs.diag<-diag(prs.s)
  
    #P1 - p1p1'
    M1<-(prs.s %*% t(prs.s))
    M2<-prs.diag - M1 
    
    #X'1(P1-p1p1')X1 = information matrix 
    X<-as.matrix(choice_set)
    M<-t(X) %*% M2  %*% X
  
    
    #sum information 
    sum_M<-sum_M + M
  
  }
 return(as.matrix(sum_M))
  }


#######################################################################################################

#Specify:-------------------------------------------------------
levels_attributes=c(3,3,3)
n_alts<-2                                                     # alternatives per set    
betas<-c(0.805, -0.080, -0.603, 0.25,-0.3, 0.70)              # Beta values             
n_sets<-10                                                    #nr of best sets to store
sum_M<-0

#---------------------------------------------------------------

#1.Generate candidate design (full design) 
candidate_profiles<-full.profiles(contrasts= "contr.sum")

k=100
D<-numeric(k)

for(i in 1:k){
#2.Generate starting design 
start<-start.design()

#3.calculate error 
INFO<-sum.information()

#D-error design
D[i]<-det(INFO)^(-1/length(betas))


}


