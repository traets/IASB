#Configure code 
rm(list=ls())

#Set working directory:
setwd("C:\Users\u0105757\Desktop\sequential_design.sas\R.code\IASB_R")

  
#Set values
nbeta<-8           #nr of betas to be estimated. 
nres<- 120      #nr of respondents 
ndraws<-512  #nr of random draws for each pp, for each beta 

#create random nr matrix 
latticefull<-as.data.frame(matrix(data=NA, nrow = nres*ndraws,ncol=nbeta))

for (i in 1:nbeta)
{latticefull[,i]<-rnorm(nres*ndraws,0,2)}

#Create beta's for each respondent 
sigma_true<- 0.25 #heterogeneity
mu_true<- as.vector(0.5*c(-1.0,0,-1.0,0,-1.0,0,-1.0,0))


beta_true<-as.data.frame(matrix(data=NA, nrow = nbeta,ncol=nres))
for (i in 1:nbeta)
{beta_true[i,]<-rnorm(nres,mu_true[i],sigma_true)}

