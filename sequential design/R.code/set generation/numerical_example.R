#Numerical example 
rm(list=ls())
#generate effect coded matrix -----------

#set values 
n_alts<-2
n_atris<-2
n_betas<-3         #sum_i (levels_atr_i - 1)

level_at1<-2
level_at2<-3

atri_1<-as.factor(seq(1,level_at1,1))
contrasts(atri_1)<-contr.sum

atri_2<-as.factor(seq(1,level_at2,1))
contrasts(atri_2)<-contr.sum



#create empty matrix 
full_design<-as.data.frame(matrix(data=NA, nrow = n_alts*n_betas, ncol=n_betas))



#Fill in values 



#Generate each possible combination
x<-c(1,2,3,4)
x<-as.factor(x)
contrasts(x)<-contr.sum
