save.image("C:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\output\\analyse_univchoice")

# regcoef labels
betalabel<-c("employability excellent job prospect vs poor",
"employability average job prospect vs poor",
"influence_parents strongly recommended by parents vs not particularly",
"influence_friends school also chosen by friends vs not",
"ranking top university ranking vs low",
"ranking average university ranking vs low",
"orientation strong academic orientation vs business",
"supervision small groups with direct contact to professors vs large")

# attribute level labels
sigmalabel<-c("employability excellent job prospect",
"employability average job prospect",
"employability poor job prospect",
"influence_parents strongly recommended by parents",
"influence_parents not particularly recommended by parents",
"influence_friends school also chosen by friends",
"influence_friends school not chosen by friends",
"ranking top university ranking",
"ranking average university ranking",
"ranking low university ranking",
"orientation strong academic orientation",
"orientation strong business orientation",
"supervision small groups with direct contact to professors",
"supervision large groups with little direct contact to professors")

library(sfsmisc)
library(abind)

logit<-function(x){log(x/(1-x))}
logist<-function(x){exp(x)/(1+exp(x))}

#############################################
# design info
############################################

## schatting multinomiaal model

P<-6 # aantal attributen
PK<-14 ## aantal attribuutlevels
J<-8 # aantal regressiecoef
M<-3 #aantal alternatieven
N<-4 #aantal alternatieven + opt-out
L<-8 # aantal mogelijke consideration sets 2^M
S<-14 # aantal keuzesets
I<-132 # aantal personen
A<-42 # aantal alternatieven =S*M
T<-1 #SC regel

# H = aantal patronen voor subsetconj regel, zie SCmat

if (T==2){
r1<-1
r2<-7}
if (T==3){
r1<-1
r2<-22}
if (T==4){
r1<-1
r2<-42}
if (T==5){
r1<-58
r2<-64}
H<-r2-r1+1

# binair design

desbin<-matrix(scan("C:\\data_HUB\\projecten\\74_university_choice_AAD\\datasets\\design_binair.txt"), nrow=42,byrow=T)

# design voor regressie met optout

desreg<-matrix(scan("C:\\data_HUB\\projecten\\74_university_choice_AAD\\datasets\\seldesign.txt"), nrow=42,byrow=T)

# choice data

choice<-read.table("C:\\data_HUB\\projecten\\74_university_choice_AAD\\datasets\\data_LGC_choice.dat",header=TRUE)
choice<-matrix(choice$choice,ncol=14,byrow=T)


#data.ISN<-array(rep(0,I*S*N),c(I,S,N))
#for (i in 1:I){
# for (s in 1:S){
#        if (choice[i,s]==1) {data.ISN[i,s,]<-c(1,0,0,0)}
#  else  if (choice[i,s]==2) {data.ISN[i,s,]<-c(0,1,0,0)}
#  else  if (choice[i,s]==3) {data.ISN[i,s,]<-c(0,0,1,0)}
#  else  if (choice[i,s]==4) {data.ISN[i,s,]<-c(0,0,0,1)}
#}}

#save(data.ISN,file="C:\\data_HUB\\projecten\\74_university_choice_AAD\\datasets\\dataISN.Rdata")
load("C:\\data_HUB\\projecten\\74_university_choice_AAD\\datasets\\dataISN.Rdata")


# patroon matrix
binmat<-matrix(rep(0,M*(2^M)),nrow=(2^M))
for (i in 0:(2^M-1)){binmat[i+1,]<-c(digitsBase(i,base=2,M))}
pat <- binmat

# SC regel
SCmat<-matrix(rep(0,P*(2^P)),nrow=(2^P))
for (i in 0:(2^P-1)){SCmat[i+1,]<-c(digitsBase(i,base=2,P))}

SCmat<-cbind(SCmat,apply(SCmat,1,sum))
ord<-order(SCmat[,P+1])
SCmat<-SCmat[ord,]

######################
# objective function
#######################


f<-function(param){

#prior
prior<-sum(log(dnorm(param)))

# bereken kans dat alternatief m in consideration set zit volgens SC-T model

if (T==1){
probc.SM<-1-matrix(exp(desbin%*%log(1-logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)}

else if (T==2|T==3|T==4){
temp<-logist(param[(J+1):(J+PK)])%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-1-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}

else if (T==5){
temp<-logist(param[(J+1):(J+PK)])%o%rep(1,S*M)
# matrix met attribuutlevelkansen per alternatief
temp2<-t(matrix(temp[t(desbin)==1],byrow=FALSE,nrow=P))

temp2.APH<-temp2%o%rep(1,H)
SCpat.APH<-rep(1,A)%o%t(SCmat[r1:r2,1:P])

#kans dat alternatief m in set s in consideration set zit
probc.SM<-matrix(apply(exp(apply(SCpat.APH*log(temp2.APH)+(1-SCpat.APH)*log(1-temp2.APH),c(1,3),sum)),1,sum),nrow=S,byrow=T)
}

else if (T==6){
probc.SM<-matrix(exp(desbin%*%log(logist(param[(J+1):(J+PK)]))),nrow=S,byrow=T)}

# bereken kans op een bepaalde consideration set per keuzeset
probcs.LS<-exp(pat%*%t(log(probc.SM))+(1-pat)%*%t(log(1-probc.SM)))

# bereken conditionele kansen
u.LSM<-rep(1,L)%o%matrix(exp(desreg%*%param[1:J]),nrow=S,byrow=T)
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

#T=1
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT1<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT1<-nlm(f,mT1$par,hessian=TRUE)
sig.T1<-as.matrix(logist(mT1$estimate[(J+1):(J+PK)]))
beta.T1<-as.matrix(mT1$estimate[1:J])
rownames(sig.T1)<-sigmalabel
rownames(beta.T1)<-betalabel
se.T1<-sqrt(diag(solve(mT1$hessian)))
dev.T1<-mT1$minimum+2*sum(log(dnorm(mT1$estimate)))
BIC.T1<-dev.T1+log(132)*length(mT1$estimate)
AIC.T1<-dev.T1+2*length(mT1$estimate)

#T=2
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT2<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT2<-nlm(f,mT2$par,hessian=TRUE)
sig.T2<-as.matrix(logist(mT2$estimate[(J+1):(J+PK)]))
beta.T2<-as.matrix(mT2$estimate[1:J])
rownames(sig.T2)<-sigmalabel
rownames(beta.T2)<-betalabel
se.T2<-sqrt(diag(solve(mT2$hessian)))
dev.T2<-mT2$minimum+2*sum(log(dnorm(mT2$estimate)))
BIC.T2<-dev.T2+log(132)*length(mT2$estimate)
AIC.T2<-dev.T2+2*length(mT2$estimate)

#T=3
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT3<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT3<-nlm(f,mT3$par,hessian=TRUE)
sig.T3<-as.matrix(logist(mT3$estimate[(J+1):(J+PK)]))
beta.T3<-as.matrix(mT3$estimate[1:J])
rownames(sig.T3)<-sigmalabel
rownames(beta.T3)<-betalabel
se.T3<-sqrt(diag(solve(mT3$hessian)))
dev.T3<-mT3$minimum+2*sum(log(dnorm(mT3$estimate)))
BIC.T3<-dev.T3+log(132)*length(mT3$estimate)
AIC.T3<-dev.T3+2*length(mT3$estimate)

#T=4
beta.s<-0.5*rnorm(J)
sigma.s<-logit(runif(PK))
param.s<-c(beta.s,sigma.s)
mT4<-optim(param.s,f,hessian=FALSE,method="BFGS")
mT4<-nlm(f,mT4$par,hessian=TRUE)
sig.T4<-as.matrix(logist(mT4$estimate[(J+1):(J+PK)]))
beta.T4<-as.matrix(mT4$estimate[1:J])
rownames(sig.T4)<-sigmalabel
rownames(beta.T4)<-betalabel
se.T4<-sqrt(diag(solve(mT4$hessian)))
dev.T4<-mT4$minimum+2*sum(log(dnorm(mT4$estimate)))
BIC.T4<-dev.T4+log(132)*length(mT4$estimate)
AIC.T4<-dev.T4+2*length(mT4$estimate)



