#############################################################################
## Generate localy optimal D-efficient designs given specifications        ##
## 1.Specify betas, nr of alts/choice set, nr attributes, levels,...       ##
## 2.Generate candidate set (full)                                         ##
## 3.Generate random start design                                          ##
## 4.Change alternatives of design to improve D.error  (mod.federov alg.)  ##
## 5.Save most efficient design                                            ##
#############################################################################

#use: 1. adjust the specify section  
#     2. select all --> run 

#clear environment
rm(list=ls())
#set.seed(1)

### Functions #########################################################################################

#Function generating candidate set (ALL profiles).
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

#Generate start design 
start.design<-function (){
  
  start<-data.frame(matrix(nrow=n_sets*n_alts, ncol=(sum(levels_attributes)-length(levels_attributes))))
  
  for(i in 1:nrow(start)){
    R<-round(runif(1, 1, nrow(candidate_profiles)))
    start[i,]<-as.numeric(candidate_profiles[R,])
  }
  
  #add setnr and alternative nr 
  set<-seq(1, n_sets, 1)
  start$set<-rep(set, each=n_alts)
  alt<-seq(1, n_alts, 1)
  start$alt<-rep(alt, n_sets)
  
  #return design
  return(start)
}

DB.design.info<- function (design)
{ 
  for (r in 1:nrow(beta_samples))
  {

  beta_r<-beta_samples[r, ]
  #design information for one beta sample
  for (c_set in 1:n_sets){
    
  
    #delete set and alt col   
    choice_set<-design[design$set==c_set,!(names(design) %in% drops)]
   
    ####information
    #create vector with all Utilities
    utilities<-as.vector(n_alts)
    for (i in 1:n_alts){utilities[i]<-sum(beta_r * choice_set[i,])}  
    
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
    
    #Save info of set
    set_info[[c_set]]<-M
    
  }
  
  #add to total 
  Bayesian_set_info<-Reduce('+', set_info)
  
}

  #mean information over betas 
  Bayesian_set_info<-Reduce('/', nrow(beta_samples))
  return(Bayesian_set_info)
  
}




#list with choice set information
set.info<- function (set) {
  
 
    
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
    info<-t(X) %*% M2  %*% X
    
    return(info)
    
  }
  
#D.error
D.error<- function (list){
sum_M<-Reduce('+', list)
D_error<-(det(sum_M))^(-1/length(betas))
return(D_error)
}


#######################################################################################################

#Specify:---------------------------------------------------------------------------------------------
levels_attributes=c(3,3,3,3)
n_alts<-2                                                               # alternatives per set    
n_sets<-15                                                              # nr of choice sets in design
n_draws<-100                                                            # nr of draws 
mu_vector<-c(0.2, 0.5, 0.06, -0.4, 1.1, -0.87, 0.1, 0.03)
var_vector<-c(1,1,1,1,1,1,1,1)

#start values:----------------------------------------------------------------------------------------
sum_M<-0                                                                 # initial value information M
drops <- c("set","alt")                                                  # add or delete col names                                                               # nr of initial starting design
converge<-F
set_info<-list()
Bayesian_set_info<-list()

#plot values
teller<-1
y<-NULL
time<-NULL
r_time<-NULL
round<-1

#libraries
library(MASS)

#-----------------------------------------------------------------------------------------------------

#1.Generate candidate design (full design) 
candidate_profiles<-full.profiles(contrasts= "contr.sum")

#3.Generate Betas 
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
beta_samples  

#2.Generate starting design 
start_des<-start.design()
work_des<-start_des


#3.calculate D-error startdesign
start_info<-DB.design.info(design=start_des) #function does not work properly yet. 

work_info<-start_info

D_start<-D.error(start_info)


#4.Change alternatives with modified federov algorithm. 
# save design that is most efficient.  

#start clock
ptm <- proc.time()

#mod.federov algorithm
while (!converge)
{
  
  #set round design
  round_des<-work_des
 
for (temp.row in 1:nrow(work_des)){

   setnr<-work_des[temp.row,]$set
  #---------------------------------------------------------------------
   for (cand in 1:nrow(candidate_profiles)){
     
     #Change alternative
     work_des[temp.row, 1:length(betas)]<-candidate_profiles[cand,]
     
     #select changed choice set
     choice_set<-work_des[work_des$set==setnr,!(names(work_des) %in% drops)]
     
     #calculate info of changed set 
     info_set<-set.info(choice_set)
     
     #change with list
     work_info[[setnr]]<-info_set
     
     #calculate error
     D<-D.error(list = work_info)
     
     #if better --> keep
      if (D < D_start)
     {
       best_row<- as.numeric(work_des[temp.row,]) 
       D_start<-D
       change<-T
       
       #for plotting
       teller<-teller+1
       y[teller]<-D
       timer<-as.list(proc.time() - ptm)
       time[teller]<-timer$elapsed
       print(D)
       flush.console()
       plot(time,y,type='l',ylab = "D.error")
       abline(v=r_time, col= "red")
      }
 
   }
  #------------------------------------------------------------------
  
  #### changes after loop through one row ##################
 
    #Keep row from temp design if improvement 
    if (change)
    {work_des[temp.row, ] <- best_row
     work_info<-design.info(work_des) }
    #Keep row from round design if no improvement  
    else
    {work_des[temp.row, ] <- round_des[temp.row,]
     work_info<-design.info(work_des) }  
                                                      
  #########################################################
  
    change<-F
    
  
  }
 #------------------------------------------------------------------
  #check if changes
  converge<-isTRUE(all.equal(work_des,round_des))
  
  #for plotting 
  #round<-round+1
  #print(round)
  
  #round_time<-as.list(proc.time() - ptm)
  #r_time[round]<-round_time$elapsed
}

#end clock
proc.time() - ptm




