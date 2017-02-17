rm(list=ls())


#libraries
library("mlogit")


#read in data
data("Heating", package = "mlogit")
 
Heating 

 H <- mlogit.data(Heating, shape="wide", choice="depvar", varying=c(3:12))
 m <- mlogit(depvar~ic+oc|0, H)
 summary(m)

s<-sum(0.071111, 0.093333, 0.636667 ,0.143333, 0.055556 )
s


mlogit.data
