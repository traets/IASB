#############################################################################
## Generate localy optimal D-efficient designs given specifications        ##
## 1.Specify betas, nr of alts/choice set, nr attributes, levels,...       ##
## 2.Generate candidate set (full)                                         ##
## 3.Generate start design                                                 ##
## 4.Change alternatives of design to improve D.error  (mod.federov alg.)  ##
## 5.Save most efficient design                                            ##
#############################################################################

#clear environment
rm(list=ls())

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
  
  start<-data.frame(matrix(nrow=n_sets*n_alts, ncol=length(betas)))
  
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

#D-error for a design
D.error<- function (design) {
  
  for ( c_set in 1:n_sets){
    
    #delete set and alt col   
    choice_set<-design[design$set==c_set,!(names(design) %in% drops)]
   
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
  
  D_error<-det(sum_M)^(-1/length(betas))
  return(D_error)
}


#######################################################################################################

#Specify:---------------------------------------------------------------------------------------------
levels_attributes=c(3,3,3,3)
n_alts<-2                                                                # alternatives per set    
betas<-c(0.25847, -0.54128, 0.15, -0.48,-0.18, 0.64, -0.12, 0.33)        # Beta values             
n_sets<-8                                                               # nr of choice sets in design
sum_M<-0                                                                 # initial value information M
drops <- c("set","alt")                                                  # add or delete col names  
n_starts<-1                                                              # nr of initial starting design
teller<-1

x<-NULL
y<-NULL
#-----------------------------------------------------------------------------------------------------

#1.Generate candidate design (full design) 
candidate_profiles<-full.profiles(contrasts= "contr.sum")

#2.Generate starting design 
start_des<-start.design()
start_start<-start_des
start<-start_des
end_design<-start_des

#3.calculate D-error startdesign
D_start<-D.error(design=start)

plot(teller, D_start)
#4.Change alternatives with modified federov algorithm. 
# save design that is most efficient.  
converge<-F

while (!converge)
{
  #teller for plot
  teller<-teller+1
  
  #end design of previous run = new start design 
  start_des<-end_design
  start<-end_design 

for (temp.row in 1:nrow(start)){

  if(temp.row>1)
  {
  #Keep row from end design if  improvement 
     if (change)
     {start[temp.row-1, ] <- end_design[temp.row-1,]}
  
  #Keep row from start design if no improvement  
     else
     {start[temp.row-1, ] <- start_des[temp.row-1,]}  
  }
  
  change<-F
  
   for (cand in 1:nrow(candidate_profiles)){
     
     #Change alternative
     start[temp.row, 1:length(betas)]<-candidate_profiles[cand,]
     
     #calculate error
     D<-D.error(design = start)
     
     #if better --> keep
     if (is.nan(D)||is.infinite(D))
     {}
     else if (D < D_start)
     {
       end_design<- as.data.frame(as.matrix(start)) 
       D_start<-D
      
       x[teller]<-teller
       y[teller]<-D
       
       change<-T
     }
    
  }}

  converge<-isTRUE(all.equal(end_design,start_des))
  plot(x,y)
}



