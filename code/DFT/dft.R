#MDFT
rm(list=ls())

#1 attribute evaluation 
Me<-c(1,2,3) #economic attribute of each alternative
Mq<-c(0.5,2,3) #quality attribute of each alternative

M<-cbind(Me,Mq)  #value matrix

#2 attention weight allocated to each attribute at particular moment t in time
weights<-c(0.3, 0.7) #probs of shifting attention 

We<-1 #prob of attending to economic attribute
Wq<-0 #prob of attending to quality attribute

W<-c(We,Wq)  #weight vector

#weighted value of each alternative at each time point 
MW<-M%*%W    #MW(t)

#3 comparison process 
C<-matrix(data=c(1, -0.5, -0.5, -0.5, 1, -0.5, -0.5, -0.5, 1), nrow= 3, ncol= 3) #contrast matrix


###valence vector 
V<-C%*%MW

### Preference: P(t+1) = SP(t) + V(t+1)

#initial preference
Po<-0

#feedback matrix
S<-matrix(data=c(0.8, -0.005, -0.005, -0.005, 0.8, -0.005, -0.005, -0.005, 0.8), nrow= 3, ncol= 3)
S2<-matrix(data=c(0.8, -0.02, -0.02, -0.02, 0.8, -0.02, -0.02, -0.02, 0.8), nrow= 3, ncol= 3)
S0<-matrix(data=c(0.8, -0, -0, -0, 0.8, -0, -0, -0, 0.8), nrow= 3, ncol= 3)

S%*%V


