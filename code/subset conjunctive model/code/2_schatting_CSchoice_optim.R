####
# Subset conjunctive model 
####


#libraries
library(sfsmisc)
library(abind)

### DATA

#levels meaning (3,3,3,4,2)
label<-c("rent 250",
	"rent 350",
	"rent 450",
	"distance 5 min",
	"distance 10 min",
	"distance 15 min",
	"noise  bad",
	"noise moderate",
	"noise very good",
	"facilities shared toilet, shower & kitchen",
	"facilities private bath & toilet, shared kitchen",
	"facilities private toilet, shared shower & kitchen",
	"facilities private toilet, shower & kitchen",
	"condition not renovated",
	"condition newly renovated")

#parameter meaning (1,2,2,3,1)
betalabel<-c("crent",
             "distance 10 min vs 5 min",
             "distance 15 min vs 5 min",
             "noise moderate vs bad",
             "noise very good vs bad",
             "facilities private bath & toilet, shared kitchen vs all shared",
             "facilities private toilet, shared shower & kitchen vs all shared",
             "facilities private toilet, shower & kitchen vs all shared",
             "condition newly renovated vs not renovated")

### FUNCTIONS
logit<-function(x){log(x/(1-x))}
logist<-function(x){exp(x)/(1+exp(x))}

############################################
# design info
############################################

## schatting multinomiaal model 
#settings
P<-5 # aantal attributen
PK<-15 ## aantal attribuutlevels
J<-9 # aantal regressiecoef
M<-5 #aantal alternatieven
N<-6 #aantal alternatieven + opt-out
L<-32 # aantal mogelijke consideration sets
S<-14 # aantal keuzesets
I<-107 # aantal personen
A<-70 # aantal alternatieven (14 choice sets * 5 alternatives)
T<-2 #SC regel

# H = aantal patronen voor subsetconj regel, zie SCmat
if (T==1){
r1<-2
r2<-32}
if (T==2){
r1<-1
r2<-6}
if (T==3){
r1<-17
r2<-32}
if (T==4){
r1<-27
r2<-32}
if (T==5){
r1<-32
r2<-32}
H<-r2-r1+1

# patroon matrix
pat<-matrix(rep(0,M*(2^M)),nrow=(2^M))
for (i in 0:(2^M-1)){pat[i+1,]<-c(digitsBase(i,base=2,M))}

# SC regel
SCmat<-matrix(rep(0,P*(2^P)),nrow=(2^P))
for (i in 0:(2^P-1)){SCmat[i+1,]<-c(digitsBase(i,base=2,P))}
SCmat<-cbind(SCmat,apply(SCmat,1,sum))
ord<-order(SCmat[,(P+1)])
SCmat<-SCmat[ord,]

# lees binair design
desbin<-matrix(scan("C:/Users/u0105757/Desktop/code/subset conjunctive model/data/design_binair.txt"),nrow=70,byrow=T)

# lees design voor regressie
# selectie van price + K_p-1 dummies per attribuut (zie betalabel)
# 3,3,3,4,2
desreg<-matrix(scan("C:/Users/u0105757/Desktop/code/subset conjunctive model/data/seldesign.txt"),nrow=70,byrow=T)

#############################
# deviance
#############################

f<-function(param){

#prior
prior<-sum(log(dnorm(param)))

# bereken kans dat alternatief m in consideration set zit volgens SC-T model

if (T==1){
probc.SM<-1-matrix(exp(desbin%*%log(1-logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)
}

else if (T==2){
temp<-logist(param[(J+1):(J+PK)])%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-1-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}


else if ((T==3)|(T==4)){
temp<-logist(param[(J+1):(J+PK)])%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}

else if (T==5){
probc.SM<-matrix(exp(desbin%*%log(logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)}

# bereken kans op een bepaalde consideration set per keuzeset
probcs.LS<-exp(pat%*%t(log(probc.SM))+(1-pat)%*%t(log(1-probc.SM)))

# bereken conditionele kansen
u.LSM<-rep(1,L)%o%matrix(exp(desreg%*%param[1:9]),nrow=S,byrow=T)
pat.LSM<-aperm(pat%o%rep(1,S),c(1,3,2))
cprob.LSM<-u.LSM*pat.LSM
temp1<-apply(cprob.LSM,c(1,2),sum)%o%rep(1,M)
temp1[1,,]<-1
cprob.LSM<-cprob.LSM/temp1
#voeg opt-out toe
coptout<-array(rep(0,L*S),c(L,S,1))
coptout[1,,]<-1
condprob.LSN<-abind(cprob.LSM,coptout,along=3)

#bereken marginale kans om een alternatief te kiezen uit een set
margprob.SN<-apply(probcs.LS%o%rep(1,N)*condprob.LSN,c(2,3),sum)

pred<-rep(1,I)%o%margprob.SN
f<--2*(sum(data.ISN*log(pred))+prior)
}

####################################
# genereer data
###################################


gendata<-function(wbeta,wsigma){

# bereken kans dat alternatief m in consideration set zit volgens SC-T model

if (T==1){
probc.SM<-1-matrix(exp(desbin%*%log(1-logist(wsigma))),nrow=S,byrow=T)
}

else if (T==2){
temp<-logist(wsigma)%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-1-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}


else if ((T==3)|(T==4)){
temp<-logist(wsigma)%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}

else if (T==5){
probc.SM<-matrix(exp(desbin%*%log(logist(wsigma))),nrow=S,byrow=T)}

# bereken kans op een bepaalde consideration set per keuzeset
probcs.LS<-exp(pat%*%t(log(probc.SM))+(1-pat)%*%t(log(1-probc.SM)))

# bereken conditionele kansen
u.LSM<-rep(1,L)%o%matrix(exp(desreg%*%wbeta),nrow=S,byrow=T)
pat.LSM<-aperm(pat%o%rep(1,S),c(1,3,2))
cprob.LSM<-u.LSM*pat.LSM
temp1<-apply(cprob.LSM,c(1,2),sum)%o%rep(1,M)
temp1[1,,]<-1
cprob.LSM<-cprob.LSM/temp1
#voeg opt-out toe
coptout<-array(rep(0,L*S),c(L,S,1))
coptout[1,,]<-1
condprob.LSN<-abind(cprob.LSM,coptout,along=3)

#bereken marginale kans om een alternatief te kiezen uit een set
margprob.SN<-apply(probcs.LS%o%rep(1,N)*condprob.LSN,c(2,3),sum)

# data genereren

data.ISN<-array(rep(0,I*S*N),c(I,S,N))
for (s in 1:S){data.ISN[,s,]<-t(rmultinom(I,1,margprob.SN[s,]))}
gendata<-list(data=data.ISN,margprob=margprob.SN)
}


# genereer parameters en data
wsigma<-logit(runif(PK))
wbeta<-0.5*rnorm(J)
genm<-gendata(wbeta,wsigma)
data.ISN<-genm$data

#####################################
# schatten model op gegenereerde data
#####################################

beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
m1<-optim(param.s,f,hessian=FALSE,method="BFGS")

par(mfrow=c(1,2))
plot(m1$par[1:9],wbeta)
plot(m1$par[10:24],wsigma)

m2<-nlm(f,m1$par,hessian=TRUE)
par(mfrow=c(1,2))
plot(m2$estimate[1:9],wbeta)
plot(m2$estimate[10:24],wsigma)



#####################################
## analyse student data
#####################################

data.ISN<-student

#T=1
T<-1
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT1<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT1<-nlm(f,mT1$par,hessian=TRUE)
sig.T1<-as.matrix(logist(mT1$estimate[10:24]))
beta.T1<-as.matrix(mT1$estimate[1:9])
rownames(sig.T1)<-label
rownames(beta.T1)<-betalabel
se.T1<-sqrt(diag(solve(mT1$hessian)))


#T=2
T<-2
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT2<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT2<-nlm(f,mT2$par,hessian=TRUE)
sig.T2<-as.matrix(logist(mT2$estimate[10:24]))
beta.T2<-as.matrix(mT2$estimate[1:9])
rownames(sig.T2)<-label
rownames(beta.T2)<-betalabel
se.T2<-sqrt(diag(solve(mT2$hessian)))

?optim

#T=3
T<-3
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT3<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT3<-nlm(f,mT3$par,hessian=TRUE)
sig.T3<-as.matrix(logist(mT3$estimate[10:24]))
beta.T3<-as.matrix(mT3$estimate[1:9])
rownames(sig.T3)<-label
rownames(beta.T3)<-betalabel
se.T3<-sqrt(diag(solve(mT3$hessian)))

#T=4
T<-4
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT4<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT4<-nlm(f,mT4$par,hessian=TRUE)
sig.T4<-as.matrix(logist(mT4$estimate[10:24]))
beta.T4<-as.matrix(mT4$estimate[1:9])
rownames(sig.T4)<-label
rownames(beta.T4)<-betalabel
se.T4<-sqrt(diag(solve(mT4$hessian)))

#T=5
T<-5
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT5<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT5<-nlm(f,mT5$par,hessian=TRUE)
sig.T5<-as.matrix(logist(mT5$estimate[10:24]))
beta.T5<-as.matrix(mT5$estimate[1:9])
rownames(sig.T5)<-label
rownames(beta.T5)<-betalabel
se.T5<-sqrt(diag(solve(mT5$hessian)))

sum(log(dnorm(mT3$estimate)))



## schatten conditional logit met opt-out

desreg2<-matrix(scan("C:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\datasets\\seldesign_optout.txt"),nrow=84,byrow=T)

choice<-c(aperm(data.ISN,c(3,2,1)))

test.data <- data.frame(
  person=gl(I,S*N,I*S*N),
  set=gl(S,N,I*S*N),
  alt=gl(N,1,I*S*N),
  choice=choice
  )

test.data <- within(test.data,{
 X1<-c(t(desreg2[,1]))
 X2<-c(t(desreg2[,2]))
 X3<-c(t(desreg2[,3]))
 X4<-c(t(desreg2[,4]))
 X5<-c(t(desreg2[,5]))
 X6<-c(t(desreg2[,6]))
 X7<-c(t(desreg2[,7]))
 X8<-c(t(desreg2[,8]))
 X9<-c(t(desreg2[,9]))
 X10<-c(t(desreg2[,10]))
})

library(mclogit)
m<-mclogit(cbind(choice,set)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=test.data)
summary(m)

############################################################################
param<-param.s
probc.SM<-1-matrix(exp(desbin%*%log(1-logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)
probc.SM<-matrix(exp(desbin%*%log(logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)

probcs.LS<-exp(pat%*%t(log(probc.SM))+(1-pat)%*%t(log(1-probc.SM)))


