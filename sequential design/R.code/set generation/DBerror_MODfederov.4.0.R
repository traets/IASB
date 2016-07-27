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
set.seed(45)

#######################################################################################################

#Specify:---------------------------------------------------------------------------------------------
levels_attributes=c(3,3,3,3)
n_alts<-2                                                               # alternatives per set    
n_sets<-15                                                              # nr of choice sets in design
n_draws<-1000                                                            # nr of draws 
mu_vector<-c(0.2, 0.5, 0.06, -0.4, 1.1, -0.87, 0.1, 0.03 )
var_vector<-c(1,1,1,1,1,1,1,1)

#start values:----------------------------------------------------------------------------------------
sum_M<-0                                                                 # initial value information M
drops <- c("set","alt")                                                  # add or delete col names                                                               # nr of initial starting design
converge<-F
design_info<-list()

#plot values-----------------------------------------------------------------------------------------
teller<-1
y<-NULL
time<-NULL
r_time<-NULL
round<-1

#libraries-------------------------------------------------------------------------------------------
library(MASS)
library(compiler)
#######################################################################################################

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

#Bayesian set information
D.error<- function (beta_r) {
 
  for (c_set in 1:n_sets)
  {
  
    #delete set and alt col   
    choice_set<-design[design$set==c_set,!(names(design) %in% drops)]
     
    ####information = X'1(P1-p1p1')X1
    Mset<-as.matrix(choice_set)
    set_info<-(t(as.matrix(choice_set)) %*% (diag(as.vector(exp(Mset%*%beta_r)/sum(exp(Mset %*% beta_r))))
             -(as.vector(exp(Mset%*%beta_r)/sum(exp(Mset %*% beta_r)))%*% t(as.vector(exp(Mset%*%beta_r)/
              sum(exp(Mset %*% beta_r)))))) %*% (as.matrix(choice_set)))
 
    #design info
    design_info[[c_set]]<-set_info
  }
  
  
  #mean information over betas 
  sum_M<-Reduce('+', design_info)
  D<-(det(sum_M))^(-1/length(beta_r))
  return(D)
}
DB.error<-mean(apply(beta_samples,1,D.error))


#mod.federov algorithm
mod.federov<-function(){while (!converge){
  
  #set round design
  round_des<-work_des
  
  for (temp.row in 1:nrow(work_des)){
    
    setnr<-work_des[temp.row, ]$set
    #---------------------------------------------------------------------
    for (cand in 1:nrow(candidate_profiles)){
      
      #Change alternative
      work_des[temp.row, 1:length(mu_vector)]<-candidate_profiles[cand,]
      
      #select changed choice set
      choice_set<-work_des[work_des$set==setnr,!(names(work_des) %in% drops)]
      
      #calculate info of changed set 
      info_set<-Bayesian.set.info2(choice_set)
      
      #change with list
      work_info[[setnr]]<-info_set
      
      #calculate error
      DB<-DB.error(design_info =  work_info)
      
      #if better --> keep
      if (DB < DB_start)
      {
        best_row<- as.numeric(work_des[temp.row,]) 
        DB_start<-DB
        change<-T
        
        #for plotting
        teller<-teller+1
        y[teller]<-DB
        timer<-as.list(proc.time() - ptm)
        time[teller]<-timer$elapsed
        print(DB)
        flush.console()
        plot(time,y,type='l',ylab = "DB.error")
        abline(v=r_time, col= "red")
      }
      
    }
    #------------------------------------------------------------------
    
    #### changes after loop through one row ##################
    
    #Keep row from temp design if improvement 
    if (change)
    {work_des[temp.row, ] <- best_row
    work_info<-Bayesian.design.info(work_des) }
    #Keep row from round design if no improvement  
    else
    {work_des[temp.row, ] <- round_des[temp.row,]
    work_info<-Bayesian.design.info(work_des) }  
    
    #########################################################
    
    change<-F
    
  }
  #------------------------------------------------------------------
  #check if changes
  converge<-isTRUE(all.equal(work_des,round_des))
  
  #plotting
  round<-round+1
  print(round)
  
  round_time<-as.list(proc.time() - ptm)
  r_time[round]<-round_time$elapsed
}
  return(round_des)
  
}

#######################################################################################################


###1.Generate candidate design (full design) 
candidate_profiles<-full.profiles(contrasts= "contr.sum")

###2.Generate Beta samples 
beta_samples<-mvrnorm(n = n_draws, mu= mu_vector, Sigma=diag(var_vector), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

###3.Generate starting design 
start_des<-start.design()
work_des<-start_des

###4.calculate DB-error startdesign
start_info<-Bayesian.design.info(design=start_des) 
work_info<-start_info
DB_start<-DB.error(start_info)

###5.Change alternatives with modified federov algorithm.
  
#start clock
ptm <- proc.time()
#ALGORITHM
DB_opt_design<-mod.federov()
#end clock
proc.time() - ptm

###6.Write design to excel file. 
write.xlsx(x=round_des, file = "C:/Users/u0105757/Desktop/sequential_design.sas/R.code/output_designs/end.design3.2.xlsx",sheetName = "TestSheet", row.names = FALSE)


