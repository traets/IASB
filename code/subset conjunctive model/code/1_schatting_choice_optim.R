## schatting multinomiaal model

P<-5
M<-3
S<-5
I<-200

## data genereren

wbeta<-0.5*rnorm(P)
x.SMP<-array(rnorm(S*M*P),c(S*M,P))
util.SM<-matrix(exp(x.SMP%*%wbeta),nrow=S,byrow=T)
kans.SM<-util.SM/apply(util.SM,1,sum)%o%rep(1,M)

data.ISM<-array(rep(0,I*S*M),c(I,S,M))
for (s in 1:S){data.ISM[,s,]<-t(rmultinom(I,1,kans.SM[s,]))}

f<-function(beta)
{
util.SM<-matrix(exp(x.SMP%*%beta),nrow=S,byrow=T)
kans.SM<-util.SM/apply(util.SM,1,sum)%o%rep(1,M)
pred<-rep(1,I)%o%kans.SM
f<--2*sum(data.ISM*log(pred))
}

beta.o<-0.5*rnorm(P)
m1<-optim(beta.o,f,hessian=TRUE,method="BFGS")


choice<-c(aperm(data.ISM,c(3,2,1)))

test.data <- data.frame(
  person=gl(I,S*M,I*S*M),
  set=gl(S,M,I*S*M),
  alt=gl(M,1,I*S*M),
  choice=choice
  )

test.data <- within(test.data,{
 X1<-c(t(x.SMP[,1]))
 X2<-c(t(x.SMP[,2]))
 X3<-c(t(x.SMP[,3]))
 X4<-c(t(x.SMP[,4]))
 X5<-c(t(x.SMP[,5])) 
})

library(mclogit)
m<-mclogit(cbind(choice,set)~X1+X2+X3+X4+X5,data=test.data)
summary(m)



