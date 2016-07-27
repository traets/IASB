rm(list=ls())

#set.seed(45)

DBS<-NULL
#Specify:#############################################################################################
levels_attributes=c(2,2,2)
n_alts<-2                                                               # alternatives per set    
n_sets<-10                                                              # nr of choice sets in design
n_draws<-1000                                                            # nr of draws 
mu_vector<-c(0.2, 0.5, 0.06)
var_vector<-c(1,1,1)
#######################################################################################################

#start values:----------------------------------------------------------------------------------------                                                                                                             # add or delete col names                                                               # nr of initial starting design
converge<-F
group<-rep(seq(1,n_sets,1),each=n_alts)  

#libraries-------------------------------------------------------------------------------------------
library(MASS)

### FUNCTIONS #########################################################################################

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
  
  return(as.matrix(full_D))
  
}

#Generate start design 
start.design<-function (){
  
  R<-round(runif((n_alts*n_sets), 1, nrow(candidate_profiles)))
  start<-candidate_profiles[R,]
  start<-as.matrix(start)
  
  #return design
  return(start)
}

#D error of design 
d.err.design<-function (beta, design){
  
  p<-t(t(work_des) * beta)
  p<-.rowSums(p,m= nrow(design),n=length(beta))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp,group), each=2)
  wp<-work_des*p
  
  info_design<-crossprod(wp,work_des)-crossprod(rowsum(wp,group))
  d_err<-det(info_design)^(-1/length(beta))
  
  return (d_err)
}

## PREPARATION ################################################################################

#generate candidate profiles
candidate_profiles<-full.profiles(contrasts= "contr.treatment")

#generate design
start_des<-start.design()
work_des<-start_des

#generate sample matrix
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector))

#calculate for all samples the D error. 
D_errors<-apply(beta_samples,1, d.err.design, design = start_des)
DB_start<-mean(D_errors)

### MOD Federov Algorithm #########################################################################################

#start clock
ptm <- proc.time()

while (!converge){
  
  #set round design
  round_des<-work_des
  
  for (temp.row in 1:nrow(work_des)){
    
    for (cand in 1:nrow(candidate_profiles)){
      
      #Change alternative
      work_des[temp.row, ]<-candidate_profiles[cand,]
      
      #calculate for all samples the D error. 
      D_errors<-apply(beta_samples,1, d.err.design , design = work_des)
      #DB error
      DB<-mean(D_errors)
      #if better --> keep
      if (DB < DB_start){
        
        best_row<- as.numeric(work_des[temp.row,]) 
        DB_start<-DB
        change<-T
      }
      
    }
    #------------------------------------------------------------------
    
    #### changes after loop through one row ##################
    
    #Keep row from temp design if improvement 
    if (change)
    {work_des[temp.row, ] <- best_row }
    #Keep row from round design if no improvement  
    else
    {work_des[temp.row, ] <- round_des[temp.row,] }  
    
    change<-F
    
  }
    #------------------------------------------------------------------
  #check if changes
  converge<-isTRUE(all.equal(work_des,round_des))
  
}
 

#end clock
proc.time() - ptm  


