#LOGIT discrete choice model-----------------------------------------------------
#following p 38-...Train
#--------------------------------------------------------------------------------

#Utility = observed component + unobserved component

#observed component= linear combination of variables 

#assumption=
#density function for each unobserved component (=random) of utility = gumbel distribution 
f<-function (Er.nj){exp(-Er.nj)*exp(-exp(-Er.nj))}

#plot distribution error component (unobserved)
Er.nj<-seq(-5,10,0.1)
dens<-f(Er.nj)
plot(Er.nj,dens)
sum(dens) #=1 

abline(v=0)


#cumulative distribution
cf<-function (Er.nj){exp(-exp(-Er.nj))}
cdens<-cf(Er.nj)
plot(Er.nj,cdens)

#variance of distribution = pi^2/6
(3.1415^2)/6 #=1.64

#by assuming this variance --> we are normalizing the scale of utility 

#the difference between two random terms with same mean = distribution with mean zero

#The difference between two extreme value variables is distributed logistic. 
#logistic distribution
e<-exp(1)
Er.nji<-seq(-5,10,0.1)
cF.log<-function (Er.nji) {e^Er.nji/(1+e^Er.nji)}
cdens.dif<-cF.log(Er.nji)
plot(Er.nji,cdens.dif)

#difference is distributed logistic
library(fExtremes)
?rgev

sam1<-rgev(n=10000, xi = 1, mu = 0, beta = 1)
sam2<-rgev(n=10000, xi = 1, mu = 0, beta = 1)


diff<-sam1-sam2

breaks<-seq(-5,10,3)

hist(diff,breaks=100, freq=FALSE)

?hist


#prob distr. extreme value (gumbel)
x=seq(-5,10,0.1)
prob<-dgev(x=seq(-5,10,0.1),mu=5,beta=1)
plot(x,prob)



#cumulative distribution
cf<-function (Er.nj){exp(-exp(-Er.nj))}
cdens<-cf(Er.nj)
plot(Er.nj,cdens)








