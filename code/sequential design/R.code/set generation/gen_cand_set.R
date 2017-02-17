

##Function generating random profiles.
#Numerical example



#set values
n_profiles<-2
n_attributes<-2
n_cand.profiles<-200

#attribute levels
levels_attributes<-as.vector(n_attributes)
levels_attributes<-c(2,3)

#create empty matrix
cand_design<-as.data.frame(matrix(data=NA, nrow = n_cand.profiles, ncol=n_attributes))

#Generate random combinations
for (r in 1:n_cand.profiles) {

  for (c in 1:n_attributes){

      x=as.factor(seq(1,levels_attributes[c],1))
      contrasts(x)<-contr.sum
      cand_design[r,c]<-sample(x=x, size=1, replace=TRUE)

}}


#Effect coding



