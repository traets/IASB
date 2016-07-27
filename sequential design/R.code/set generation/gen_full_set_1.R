####################################
## Initial static stage
## 1. generate All profiles
## 2. Select type of coding
####################################

#clear environment
rm(list=ls())

##Function generating ALL profiles.----------------------
rand.profiles<-function (levels_attributes) {

  #as list
  levels_attributes<-as.list(levels_attributes)

  #Make list with all atribute levels
  fun<-function (x){return(1:x)}
  levels_attributes<-lapply(X = levels_attributes, fun)

  #generate every combination
  full_design<-as.data.frame(expand.grid(levels_attributes))

  #create factors:
  col_names <- names(full_design)
  full_design[,col_names] <- lapply(full_design[,col_names] , factor)

  #generate design according to type of contrast
  list.attris<-list(colnames(full_design))
  contr.fun<-function(x){paste(x,= "contr.sum")}
  contrast.list<-lapply(X = list.attris, FUN = contr.fun)
  full_design<-as.data.frame(model.matrix(~ ., full_design,contrasts = contrast.list))

  return(full_design)

  }

full_design<-rand.profiles(levels_attributes=c(2,4,5,3))
contrast.list
col_names
list.attris
