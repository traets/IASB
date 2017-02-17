
#libraries
library(mlogit)

#Data
data("Heating", package = "mlogit")

#Run model
H <- mlogit.data(Heating, shape="wide", choice="depvar", varying=c(3:12))
m <- mlogit(depvar~ic+oc|0, H)
summary(m)


#How closely do the average probabilities match the shares
#of customers choosing each alternative ?
apply(fitted(m, outcome=FALSE), 2, mean)

#Ratio of operating cost coeficcient to the installation cost coefficient.
#--> willingness to pay. 
coef(m)["oc"]/coef(m)["ic"]

#put a contraint on model: lcc= ic + oc/0.12. 
#estimate model with this var.

H$lcc=H$ic+H$oc/0.12
mlcc <- mlogit(depvar~lcc|0, H)
mlcc
lrtest(m, mlcc)

#test significance
qchisq(0.05, df = 1, lower.tail = FALSE)

#Add alternative-specific constants to the model (leaving out"|0")
mc <- mlogit(depvar~ic+oc, H, reflevel = 'hp')
summary(mc)

#How closely do the average probabilities match the shares
#of customers choosing each alternative ?
apply(fitted(mc, outcome=FALSE), 2, mean) #exactly 

#Enter sociodemographic variables. 
#installation cost dived by income
mi <- mlogit(depvar~oc+I(ic/income), H)
summary(mi)
#--> worse

mi2 <- mlogit(depvar~oc+ic|income, H, reflevel="hp")
summary(mi2)

lrtest(mc,mi2)

#Try own model------------------------------------------------------
#mlogit(Depvar ~ generic | altspec)

mod<-mlogit(depvar~oc+ic|0, H)
mod2<-mlogit(depvar~oc+ic|rooms+region, H)
summary(mod)
summary(mod2)

#fit with frequencys 
apply(fitted(mod2, outcome=FALSE), 2, mean) 


#Use model for prediction----------------
X <- model.matrix(mc)
head(X)
alt <- index(H)$alt
chid <- index(H)$chid
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- tapply(eXb, chid, sum)
P <- eXb / SeXb[chid]
P <- matrix(P, ncol = 5, byrow = TRUE)
head(P)
apply(P, 2, mean)


#-------------------------------------------------------------------------------------
#Mixed Logit__________________________________________________________________________
#-------------------------------------------------------------------------------------
rm(list=ls())

#Libraries
library("mlogit")

#Data
data("Electricity", package = "mlogit")

#reshape
Electr <- mlogit.data(Electricity, id="id", choice="choice",
                          varying=3:26, shape="wide", sep="")


Elec.mxl <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr,
                       rpar=c(pf='n', cl='n', loc='n', wk='n', tod='n', seas='n'),
                       R=100, halton=NA, print.level=0, panel=TRUE)
summary(Elec.mxl)



