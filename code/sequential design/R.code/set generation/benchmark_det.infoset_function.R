rm(list=ls())
set.seed(45)


#Specify:#############################################################################################
levels_attributes=c(3,3,3,3)
n_alts<-2                                                               # alternatives per set    
n_sets<-15                                                              # nr of choice sets in design
n_draws<-10                                                        # nr of draws 
mu_vector<-c(0.2, 0.5, 0.06, -0.4, 1.1, -0.87, 0.1, 0.03 )
var_vector<-c(1,1,1,1,1,1,1,1)
#######################################################################################################

#start values:----------------------------------------------------------------------------------------
sum_M<-0                                                                 # initial value information M
drops <- c("set")                                                  # add or delete col names                                                               # nr of initial starting design
converge<-F

#plot values-----------------------------------------------------------------------------------------
teller<-1
y<-NULL
time<-NULL
r_time<-NULL
round<-1

#libraries-------------------------------------------------------------------------------------------
library(MASS)
library(compiler)
library(benchmark)
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
info.set<-function (choice_set, beta){
  
  x<-as.matrix(choice_set)
  x
  p <- c(exp(x %*% beta))
  p <- p / sum(p)
  set_info<-(t(x) %*% (diag(p) - tcrossprod(p)) %*% (x))
  return(set_info)
}
info.set_c <- compiler::cmpfun(info.set)


#D.err 
D.err<- function (beta, work_des) {
  
  design_info<-by(subset(work_des, select=-c(set)), work_des$set, info.set_c, beta = beta)
  design_info<-Reduce('+', design_info)
  D_err<-((prod(diag(chol(design_info))))^2)^(-1/ncol(beta_samples))
  
  return(D_err)
}
D.err_c <- compiler::cmpfun(D.err)


D.err2<- function (beta, work_des) {
  
  design_info<-by(subset(work_des, select=-c(set)), work_des$set, info.set_c, beta = beta)
  design_info<-Reduce('+', design_info)
  D_err<-det(design_info)^(-1/ncol(beta_samples))
  
  return(D_err)
}
D.err2_c <- compiler::cmpfun(D.err2)


#generate candidate profiles
candidate_profiles<-full.profiles(contrasts= "contr.treatment")

#generate design
start_des<-start.design()
work_des<-start_des

#generate sample matrix
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector))


#NEw FUNCTION
d.err.design<-function (beta, work_des){
  
  group<-rep(seq(1,n_sets,1),each=n_alts)  
  
  x<-as.matrix(subset(work_des, select=-c(set)))
  
  p<-x %*% diag(beta)
  p<-.rowSums(p,m= n_sets*n_alts,n=length(beta))
  p_total<-rep(rowsum(exp(p),group), each=2)
  p<-exp(p)/p_total
  
  
  
  info_design<-crossprod(x*p,x)-crossprod(rowsum(x*p,group))
  d_err<-det(info_design)^(-1/length(beta))
  return (d_err)
}

d.err.design_c <- compiler::cmpfun(d.err.design)

library(microbenchmark)

#calculate for all samples the D error. 
beta<-c(3,2,1,0,1,2,3,2)


microbenchmark(DB_start<-D.err(beta, work_des ),times = 1000)

microbenchmark(DB_start2<-D.err_c(beta, work_des),times = 1000)

microbenchmark(DB_start3<-D.err2(beta, work_des),times = 1000)

microbenchmark(DB_start4<-D.err2_c(beta, work_des), times = 1000)

microbenchmark(DB_start5<-d.err.design(beta, work_des), times = 1000)

microbenchmark(DB_start6<-d.err.design_c(beta, work_des), times = 1000)






