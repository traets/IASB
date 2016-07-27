rm(list=ls())

#Specify:#############################################################################################
levels_attributes=c(3,3,3,3)
n_alts<-2                                                               # alternatives per set    
n_sets<-15                                                              # nr of choice sets in design
n_draws<-100                                                          # nr of draws 
mu_vector<-c(0.2, 0.5, 0.06, -0.4, 1.1, -0.87, 0.1, 0.03 )
var_vector<-c(1,1,1,1,1,1,1,1)
#######################################################################################################

### Functions #########################################################################################

#Function generating candidate set (ALL profiles).
full.profiles<-function (contrasts) {
  
  
  #Make list with all atribute levels
  levels_attributes<-lapply(X = as.list(levels_attributes), function (x) (1:x))
  
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

#Generate start design 
start.design<-function (){
  
  R<-round(runif((n_alts*n_sets), 1, nrow(candidate_profiles)))
  start<-candidate_profiles[R,]
  
  #add setnr 
  set<-seq(1, n_sets, 1)
  start$set<-rep(set, each=n_alts)
  
  #return design
  return(start)
}

#info choice set 
info.set<-function (choice_set){
  
  x<-as.matrix(choice_set)
  set_info<-(t(x) %*% (diag(as.vector(exp(x%*%betas)/sum(exp(x %*% betas))))
                       -(as.vector(exp(x%*% betas)/sum(exp(x %*% betas)))%*% t(as.vector(exp(x%*%betas)/
                                                                                           sum(exp(x %*% betas)))))) %*% (x))
  return(set_info)
}

#D.err 
D.err<- function (work_des) {
  
  design_info<-by(subset(work_des, select=-c(set)), work_des$set, info.set)
  design_info<-Reduce('+', design_info)
  D_err<-det(design_info)^(-1/ncol(beta_samples))
  
  return(D_err)
  }

### Algorithm #########################################################################################

#generate candidate profiles
candidate_profiles<-full.profiles(contrasts= "contr.sum")

#generate design
work_des<-start.design()

#generate sample matrix
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector))

#vector to store D errors
D_errors<- numeric(n_draws)

#calculate for all samples the D error. 
for (i in 1:nrow(beta_samples)){
betas<-beta_samples[i,]
D_errors[i]<-D.err(work_des)
}

D_errors

#alternative 
test<-apply(beta_samples,1,D.err)


