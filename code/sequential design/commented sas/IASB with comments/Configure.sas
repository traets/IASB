/*libname is used for assigning the directory "C:\..." as a SAS library. This directory will have all your inputs and you will
store all you outputs -the ones you require- to this directory.
This assigns a temporary name "seqdb" to that directory. To call an input or to store an output, you should use "seqdb.input"
or "seqdb.output". You can use any name instead of "seqdb"*/
libname seqdb 'C:\Users\u0105757\Desktop\sequential design\commented sas\IASB with comments';

%let df = 4;/*number of betas to be estimated. For 4 attributes,you have to estimate 8 betas for effect-coded design matrix*/
%let nres = 1;/* number of respondents*/
/*These parameters are used for generating lattice points for random draws
this setting will give you 2^9=512 random lattice points
These lattice points are used for generating normally distributed random draws
These random draws are necessary for approximating integrals in Information matrices and/or for generating Bayesian prior parameters */
%let a = 1571;
%let b = 2;
%let m = 8; /*This will generate 2^8 = 256 random draws*/
/*Generating random lattice points*/
proc iml;

do res=1 to &nres;

start base(num);
a1=j(&m,1,0);
c1=j(&m,1,0);
a1[&m]=mod(num,&b);
c1[&m]=a1[&m];
do j=&m-1 to 1 by -1;
tem=(num-c1[j+1])/(&b**(&m-j));
a1[j,1]=mod(tem,&b);
c1[j]=c1[j+1]+a1[j]*(&b**(&m-j));
end;
return(a1);
finish;
n=&b**&m;
u=uniform(repeat(0,1,&df));
av=j(&df,1,1);
do i=2 to &df;
av[i]=mod((av[i-1]*&a),n);
end;
e=j(n,&df,0);
seq=j(&m,1,0);
kk=-&m;
do k=1 to &m;
seq[k]=&b**kk;
kk=kk+1;
end;
do i=1 to n;
ei=(seq`*base(i-1))*av`+u;
e[i,]=ei-int(ei);
end;
latt=j(n,&df,0);
do i=1 to &df;
do j=1 to n;
latt[j,i]=1-abs(2*e[j,i]-1);
end;
end;
lattice=probit(latt);

latticefull = latticefull // lattice;
end;

create seqdb.latfull from latticefull; append from latticefull;
quit;
/*The following parameters are used in D-error computation and used for simulating choice data */
libname seqdb 'C:\Users\u0105757\Desktop\sequential design\commented sas\IASB with comments';

%let df = 4;
%let nres = 1;

proc iml;

mutrue = 0.5*{-1.0,0,-1.0,0};

create seqdb.mutrue from mutrue; append from mutrue;

sigmatrue = 1.0*I(&df);
create seqdb.sigmatrue from sigmatrue; append from sigmatrue;

chol=t(root(sigmatrue));
rnorm = j(&df,&nres,.);
call randgen(rnorm,'normal');
betatrue = mutrue+(chol*rnorm);

print betatrue;
print seqdb;
create seqdb.betatrue from betatrue; append from betatrue;

quit;

%let df = 4;
/*These parameters are used in the estimation procedure as initial (starting) values*/
proc iml;

pmean=j(&df,1,0);
pcov=&df*I(&df);

create seqdb.pmean from pmean; append from pmean;
create seqdb.pcov from pcov; append from pcov;

quit;
