logit<-function(x){log(x/(1-x))}

wparam<-scan("C:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\X_test\\wparam.txt")

samp<-matrix(scan("C:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\X_test\\sample.txt"),nrow=50000,byrow=T)

par(mfrow=c(3,3))
for (i in 1:9){plot(samp[,i],ylim=c(-3,3))}

par(mfrow=c(4,4))
for (i in 10:24){plot(samp[,i],ylim=c(0,1))}


msamp<-apply(samp,2,mean)
vsamp<-sqrt(apply(samp,2,var))

par(pty="s")
par(mfrow=c(1,2))
plot(wparam[1:9],msamp[1:9])
plot(logit(wparam[10:24]),logit(msamp[10:24]))

cor(wparam[1:9],msamp[1:9])
cor(logit(wparam[10:24]),logit(msamp[10:24]))

logist<-function(x){exp(x)/(1+exp(x))}

margprob<-matrix(scan("c:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\X_test\\margprobsm.txt"),nrow=14,byrow=T)

gendat<-matrix(scan("c:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\X_test\\gendat.txt"),nrow=100,byrow=T)

temp<-t(apply(gendat[,2:15],2,table))/5000

plot(temp,margprob)


