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
  D_err<-det(design_info)^(-1/ncol(beta_samples))
  
  return(D_err)
  }



#generate candidate profiles
candidate_profiles<-full.profiles(contrasts= "contr.treatment")

#generate design
start_des<-start.design()
work_des<-start_des

#generate sample matrix
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector))

#calculate for all samples the D error. 
D_errors<-apply(beta_samples,1, D.err, work_des = work_des)
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
      work_des[temp.row, 1:length(mu_vector)]<-candidate_profiles[cand,]
      
      #calculate for all samples the D error. 
      D_errors<-apply(beta_samples,1, D.err, work_des = work_des)
      #DB error
      DB<-mean(D_errors)
      #if better --> keep
      if (DB < DB_start)
      {
        best_row<- as.numeric(work_des[temp.row,]) 
        DB_start<-DB
        change<-T
#         
#         #for plotting
#         teller<-teller+1
#         y[teller]<-DB
#         timer<-as.list(proc.time() - ptm)
#         time[teller]<-timer$elapsed
          print(DB)
#         flush.console()
#         plot(time,y,type='l',ylab = "DB.error")
#         abline(v=r_time, col= "red")
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
    
    #########################################################
    
    change<-F
    
  }
  #------------------------------------------------------------------
  #check if changes
  converge<-isTRUE(all.equal(work_des,round_des))
  
#   #plotting
#   round<-round+1
#   print(round)
#   
#   round_time<-as.list(proc.time() - ptm)
#   r_time[round]<-round_time$elapsed
}
  
#end clock
proc.time() - ptm  








