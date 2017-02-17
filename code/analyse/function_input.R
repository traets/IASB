rm(list = ls())
require(maxLik)
require(data.table)
require(sandwich)

 
detachAllData <-
function () 
{
    pos.to.detach <- (1:length(search()))[substring(search(), 
        first = 1, last = 8) != "package:" & search() != ".GlobalEnv" & 
        search() != "Autoloads" & search() != "CheckExEnv" & search() != "tools:rstudio" & search() != "TempEnv"]
    for (i in 1:length(pos.to.detach)) {
        if (length(pos.to.detach) > 0) {
            detach(pos = pos.to.detach[1])
            pos.to.detach <- (1:length(search()))[substring(search(), 
                first = 1, last = 8) != "package:" & search() != 
                ".GlobalEnv" & search() != "Autoloads" & search() != 
                "CheckExEnv" & search() != "tools:rstudio" & 
                search() != "TempEnv"]
        }
    }
}

detachAllData()

modeloutput=function(model)
{
 est<<-model$estimate
 varcov=vcov(model)
 meat1=meat(model)
 bread1=bread(model)
 meat1[is.na(meat1)]=0
 bread1[is.na(bread1)]=0
 robvarcov=sandwich(model,bread1,meat1)
 se=sqrt(diag(varcov))
 robse=sqrt(diag(robvarcov))
 trat_0=est/se
 robtrat_0=est/robse
 trat_1=(est-1)/se
 robtrat_1=(est-1)/robse
 se[model$fixed]=NA
 robse[model$fixed]=NA
 trat_0[model$fixed]=NA
 robtrat_0[model$fixed]=NA
 trat_1[model$fixed]=NA
 robtrat_1[model$fixed]=NA
 varcov[model$fixed,]=NA
 varcov[,model$fixed]=NA
 robvarcov[model$fixed,]=NA
 robvarcov[,model$fixed]=NA
 finalLL=model$maxim
 iterations=model$iterations
 zeroLL=sum(loglike(0*beta))
 initLL=sum(loglike(beta))
 params=length(beta)-sum(model$fixed)
 rho2zero=1-finalLL/zeroLL
 adjrho2zero=1-(finalLL-params)/zeroLL
 est=round(est,4)
 se=round(se,4)
 trat_0=round(trat_0,2)
 trat_1=round(trat_1,2)
 robse=round(robse,4)
 robtrat_0=round(robtrat_0,2)
 robtrat_1=round(robtrat_1,2)
 output=t(rbind(est,se,trat_0,trat_1,robse,robtrat_0,robtrat_1))
 
 cat("Model diagnosis:",model$message,"\n\n")

 cat("LL: ",finalLL,"\n\n")
 
 cat("Estimates:\n")
 print(output)

 varcov=signif(varcov,4)
 robvarcov=signif(robvarcov,4)
 
 cat("\n\nCovariance matrix:\n")
 print(varcov)

 cat("\n\nRobust covariance matrix:\n")
 print(robvarcov)
}

runmodel=function()
{
  model=maxLik(loglike,start=beta,fixed=fixedparams,method="BFGS",print.level=3)
  modeloutput(model)
}

shuffle=function(inv){
out=inv[rank(runif(length(inv)))];
out}

mlhs=function(N,d,i){
temp=seq(0,N-1)/N;
out=matrix(0,N*i,d);
j=1;
k=1;
while(j<i+1){
k=1;
while(k<d+1){
out[(1+N*(j-1)):(N*j),k]=shuffle(temp+runif(1)/N);
k=k+1}
j=j+1}
out}

halton=function(n,d){
prime=2;
out=(haltonsequence(prime,n));
i=2;
while(i<d+1){
k=0
while(k<1){
prime=prime+1;
if(sum(prime/1:prime==prime%/%1:prime)==2) k=1;
}
out=cbind(out,haltonsequence(prime,n));
i=i+1}
out}

haltonelement=function(prime,element){
H=0;
power=(1/prime);
while(element>0){
digit=(element%%prime);
H=H+digit*power;
element=element%/%prime;
power=power/prime}
H}

haltonsequence=function(prime,lengthvec){
i=1;
out=0;
while(i<lengthvec+1){
out[i]=haltonelement(prime,i);
i=i+1}
out}
