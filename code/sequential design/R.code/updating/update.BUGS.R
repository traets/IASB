###################################
# Concepts in Bayesian Inference  #
# Assignment                      #
###################################


# set working directory 
rm(list=ls())
setwd("C:/Users/u0105757/Desktop/london course/R practicals/MNL practical")

#libraries
library(xlsx)
library(R2OpenBUGS)
library(MASS) 
library(coda)
library(BRugs)
library(HDInterval)
library(dplyr)

#load DATA
df<-read.table("danish_data.dat",header=T)
df<-select(df[1:50,], choice: tc2)
df$choice<-df$choice-1

#data to openbugs
data<-list("choice"=df$choice,
           "tt1"=df$tt1,
           "tt2"=df$tt2,
           "tc1"=df$tc1,
           "tc2"=df$tc2,
           "N"=nrow(df))


#Specify parameters with initial starting values     
params <- c("interc", "beta1", "beta2")

# Starting values 
inits <-NULL

#specify model
model<-function() {
  
  #likelihood
  for (i in 1:N) {
    
    
    choice[i] ~ dbern(mu[i])
    mu[i]<-U2/U1+U2
    # define utility functions
    U1<-exp(interc+beta1*tc1+beta2*tt1)
    U2<-exp(interc+beta1*tc2+beta2*tt2)
    
  }
  
   #vague priors
   interc ~ dnorm(0, 1.0E-6)
   beta1 ~ dnorm(0, 1.0E-6)
   beta2 ~ dnorm(0, 1.0E-6)
  
}


## RUN MODEL #######################################
#Write model to Winbugs model file
model.file <- file.path(tempdir(),"model.txt") 
write.model(model, model.file)

#Run Bugs (change inits to inits2 for other initial values)
out <- bugs(data, inits, params, model.file, codaPkg=TRUE, n.iter=10000,n.thin = 1, n.chains = 2, debug= T)
#save results
out.coda <- read.bugs(out) 
#####################################################

## Evaluate convergence #####
# trace plots + posteriors 
plot(out.coda)
# autocorrelation 
acfplot(out.coda)
# gelman diagnostic
gelman.diag(out.coda)
gelman.plot(out.coda)


## summary ####################
#mean, mode and sd        
summary(out.coda)

#HPD intervals 
HPDinterval(out.coda, prob = 0.95)
hdi(out.coda, credMass = 0.95)


#### frequentist 
# MLE logistic regression:
fit.freq = glm( TSH.CAT ~ GENDER * BTHWGT,data=df, family=binomial(logit) )
show(fit.freq )

