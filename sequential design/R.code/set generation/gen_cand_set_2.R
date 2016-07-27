####################################
## Initial static stage
## 1. generate random profiles
## 2. change to effect coding
####################################

#clear environment
rm(list=ls())

##Function generating random profiles.----------------------
rand.profiles<-function (n_attributes, n_cand.profiles, levels_attributes, effectcoding) {

#create empty matrix
cand_design<-as.data.frame(matrix(data=NA, nrow = n_cand.profiles, ncol=n_attributes, effectcoding))

if (effectcoding){
betas<-sum(levels_attributes)-(length(levels_attributes))
cand_eff.design<-as.data.frame(matrix(data=NA, nrow = n_cand.profiles, ncol=betas))}

#Generate random combinations
for (r in 1:n_cand.profiles) {

  if(effectcoding){

    for (b in 1:n_attributes){

    contrasts(x)<-contr.sum
    x.con<-contrasts(x)
    cand_eff.design[r,]<-x.con[sample(x=x, size=1, replace=TRUE), ])
    }

  }else{

  for (c in 1:n_attributes){

      x=as.factor(seq(1,levels_attributes[c],1))
      cand_design[r,c]<-sample(x=x, size=1, replace=TRUE)
  }
    return(cand_design)

 }


}


sample(x=x[,1], size=1, replace=TRUE)


#create example
design<-rand.profiles(n_attributes = 4, n_cand.profiles = 200, levels_attributes = c(2,3,5,4))


#Effect coding---------------------------------------------

betas<-sum(levels_attributes)-(length(levels_attributes))
effect.design<-as.data.frame(matrix(data=NA, nrow = n_cand.profiles, ncol=betas))


x<-c(1,2,3,4)
x<-as.factor(x)
contrasts(design)<-contr.sum


#effect codes
codes<-array(0, dim=c(n_attributes,max(levels_attributes),betas))


n_attributes=2
levels_attributes = c(2,3)

for( a in 1:n_attributes){


  f<-as.factor(seq(1,levels_attributes[a],1))
  contrasts(f)<-contr.sum
  codes[a, , ]<-f

  }



for (j in 1:n_cand.profiles){

  for (k in 1:betas)  {



  }
}


x<-c(1,2,3,4)
x<-as.factor(x)
contrasts(x)<-contr.sum

x=as.factor(seq(1,levels_attributes[c],1))
contrasts(x)<-contr.sum

x.con<-contrasts(x)
x.con


x.con<-as.matrix(x)
x.con
x[1,]
cand_design[r,c]<-sample(x=x, size=1, replace=TRUE)


x[2,]




