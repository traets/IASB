# read data and make call to run for screening rules
#  note you must change path in the file names below to point to the correct
#  location

rm(list=ls())

source("rScreen.R",local=TRUE)
dyn.load("screen.dll", mode="wb")


y=read.table("camchoice.txt",header=FALSE) # the file "camchoice.txt" contains the dependent variable (y=1 or 0)
x=read.table("camdesign.txt",header=FALSE) # the file "camdesign.txt" contains the design matrix (dummy variable coding)
inp2=read.table("camatt3.txt",header=FALSE) # the file "camatt3.dat" contains the design matrix (levels coding)
ind <- matrix(scan("tindex.txt",0), ncol=3, byrow=TRUE) # the file "tindex.txt" contains the theta index for each attribute

nhh=302 # nhh = number of respondents (households)
nset=14 # nset = number of choice tasks
nsize=7 # nsize = number of alternatives per choice task
nxvar=18 # nxvar = number of attribute levels, dim(beta)
natvar=11 # natvar = number of attributes
ntheta=36 # ntheta = total number of grid points needed for discrete attributes

y = array(t(y), dim=c(nsize,nset,nhh))
X = array(t(x), dim=c(nxvar,nsize,nset,nhh))
xatt=array(t(inp2), dim=c(natvar,nsize,nset,nhh))

Data=list(y=y,X=X,xatt=xatt,ind=ind,nhh=nhh,nset=nset,
          nsize=nsize,nxvar=nxvar,natvar=natvar,ntheta=ntheta)

nu=nxvar+5
Prior=list(nu=nu,V0=nu*diag(rep(1,nxvar)),
      betabarbar=as.vector(rep(0,nxvar)),
      Abeta=.01*diag(rep(1,nxvar)))
Mcmc=list(R=1000,keep=1)

out=rScreen(Data,Prior,Mcmc)

rScreen

