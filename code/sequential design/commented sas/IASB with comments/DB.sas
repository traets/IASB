libname seqdb 'C:\Users\u0105757\Desktop\sequential design\commented sas\IASB with comments';

%let nalts = 2; /* # of alternatives*/
%let niset = 5; /* # of initial choice sets*/
%let tset = 6; /* # of total choice sets*/

%let nres = 1; /* # of respondents; of course 10 respondents won't be enough tohave accurate 
				estimates; mostly it is chosen more than 50; as 200 or 250*/
%let df = 4; /* # of betas*/
%let dfv = 8; /* parameter of the multivariate student t-distrinbution; used in Importance
				Sampling */

%let ndraws = 16;/* # of random parameter draws. For lattice points, this number is always 
				 set to 2^m, here m=4, mostly it is used as m=8 or 9 */


/*This following part generates the full-factorial design for the specified design setting
For a small design set-up, 3^2/2/10 is used. 2 attributes with 2 levels, 2 alternatives and 
10 choice sets, 5 of them generated as initial sets.*/
proc plan ordered;
factors x1=3 x2=3 /noprint;
output out=candidate;
run;
/*This gives the effect-coded version of the full-factorial design*/
proc transreg design data=candidate;
model class(x1 x2 /effects);
output out=candid;
run;

/*SAS/IML is a matrix language, starts with "proc iml" and ends with "quit"*/
proc iml WORKSIZE=6000000;
options nonotes;

use seqdb.initialsets_DB; read all into initialsets;
print initialsets;

use seqdb.betatrue; read all into truebeta; 
use seqdb.mutrue; read all into pmean; /*prior mean*/
use seqdb.sigmatrue; read all into pcov; /*prior covar*/
invcovar = inv(pcov);

use candid(keep = &_trgind); read all into candidat;
/*use candidate; read all into candidat;*/
do i = 1 to nrow(candidat)-1;
	do j = 1 to nrow(candidat)-i;
    	holddesignfull = holddesignfull // candidat[i, ];
    	holddesignfull = holddesignfull // candidat[i + j, ];
	end;
end;

/*This function is for simulating the choice data - y-values*/
start choice(beta,Xlast,y1,cset,y);
y=j(cset,1,0);
if cset=1 then y=y;
else do;
y[1:cset-1,1]=y1;
end;
p=j(&nalts,1,0);
ran=uniform(0);
utis=exp(Xlast*beta);
p=utis/sum(utis);
  if ran<=p[1] then y[cset,1]=1;
  else if ran<=p[1]+p[2] then y[cset,1]=2;
  else y[cset,1]=3;
finish;

/*This function is for generating random draws from the "Multivariate student t-distribution" 
which is used in Importance Sampling
See link for probability density function of Multivariate student t-distribution:
https://en.wikipedia.org/wiki/Multivariate_t-distribution*/
start RANDMVT(N, DF, Mean, Cov, RANZ);	
	/* check parameters */
	if N<1 then do;
	print "The requested number of observations should be at least 1:" N; stop;
	end; 
	if DF<1 then do;
	print "The degrees of freedome should be at least 1:" DF; stop;
	end;
	mMean = colvec(Mean);
	p = nrow(mMean); 
	X = j(N,p,.);
	*Z = j(p,1,.);
	A = root( Cov );
	inv = 0; 
	do i=1 to N;
		call randgen(inv,'GAMMA',DF/2);   /* Generate inv~ Gamma(DF/2,1)*/
		invW = inv/(DF/2);      /* invW ~ Gamma(DF/2, 2/DF ) */
	 	W = 1/invW;    /* W ~ Inverse Gamma(DF/2, DF/2) */
		*call randgen( Z, 'NORMAL'); 
		Z=RANZ[,i];
		r = mMean + sqrt(W)* A` * Z ; /*A` Z ~ the standard univariate normal(if p = 1) or Multivariate normal(if p>1) */
		X[i, ] = r`;
	end;
	return(X);
finish;

/*The following functions (prob, prob1, prob2, infologit and infologit2 are used in Importance
Sampling*/
/*This function is for calculating the probabilities of alternatives based on MNL model
for all design*/
start prob(X,b,cset);
p=j(&nalts,cset,0);
	do ip=1 to cset;
	xx=X[(ip-1)*&nalts+1:ip*&nalts,];
	m = max(xx*b);
	u=(xx*b)-m;
  	psub=exp(u)/sum(exp(u));
   	p[,ip]=psub;
	end;
return(p);
finish;

/*This function is for calculating the probabilities of alternatives based on MNL model
for one choice set*/
start prob1(X,b);
m = max(X*b);
u=(X*b)-m;
  psub=exp(u)/sum(exp(u));
return(psub);
finish;

/*This function is for calculating the probabilities of selected alternatives for all design;
in other words, it computes the likelihood of MNL/CL model*/
start prob2(yn,X,b);
row=nrow(X);
cset=row/&nalts;
p=1;
	do ip=1 to cset;
	xx=X[(ip-1)*&nalts+1:ip*&nalts,];
	m = max(xx*b);
	u=(xx*b)-m;
  	psub=exp(u)/sum(exp(u));
  	p=p*psub[yn[ip]];
	end;
return(p);
finish;

/*This funcion is for computing the Fisher Information Matrix of MNL model for all design*/
start infologit(X,b,cset,info1);
info1=j(&df,&df,0);
p=prob(X,b,cset);
do in=1 to cset;
xx=X[(in-1)*&nalts+1:in*&nalts,];
pp=diag(p[,in]);
   info1=info1+xx`*(pp-p[,in]*p[,in]`)*xx;
end;
finish;

/*This funcion is for computing the Fisher Information Matrix of MNL model for a choice set*/
start infologit2(X,b,info2);
p1=prob1(X,b);
   pp=diag(p1);
info2=X`*(pp-p1*p1`)*X;
finish;

/*This function is for computing the probability density function of "Multivariate Student 
t-distribution" which is used in Importance Sampling*/
start MVST(n,beta,meanbeta,cov);
dif=beta-meanbeta;
invcov=inv(cov);
diff=t(dif)*invcov*dif;
iMVSTd=1/(det(cov)**(0.5))*(1+((1/n)*diff))**(-(n+&df)/2);
return(iMVSTd);
finish MVST;

/************************main loop*******************************/

obs=&nres#&tset#&nalts;
/*Xmat and ymar are storage matrices for the sequential designs to be generated*/
Xmat=j(obs,&df,0);
ymat=j(&nres,&tset,0);

use seqdb.latfull; read all into zfull;/*Random parameter draws*/

/*this the loglikelihood of the priors, note that we assume betas are normally distributed.*/
logprior1=-&df/2*log(2*3.1415926)-0.5*log(det(pcov));
prior1=1/((2*3.1415926)**(&df/2)*det(pcov)**(0.5));

do res=1 to &nres; /*for-loop of each respondents starts here*/
holddesign = holddesignfull;
nrowcand=nrow(holddesign);
beta=truebeta[,res];
points = (res-1)*&ndraws+1:res*&ndraws; 
z = zfull[points,];
z=z`;
/*--------------------------------------------------------------*/
/*This following part is only used when there is no initial set.*/
/*(1)*/
/*
chol=t(root(pcov)); 
priordraws=pmean+chol*z;  
derrpriorsmall=10000;
nsetcand=nrowcand/&nalts;
do u=1 to nsetcand;
	l=(u-1)*&nalts+1:u*&nalts;
	X=holddesign[l,];
	derrprior = 0;
	do d = 1 to &ndraws;
	b = priordraws[,d];
	call infologit2(X,b,info);
	binfo = info + invcovar;
	detinfo = det(binfo);
	derr = detinfo##(-1/&df);
	derrprior = derrprior + derr;
	end;
	derrprior = derrprior/&ndraws;
if derrprior < derrpriorsmall then;
do;
derrpriorsmall = derrprior;
startset=X;
end;
end; 
Xorg=startset;*/
/*--------------------------------------------------------------*/

nruns=&nalts*&niset;
ss=(res-1)*nruns+1:res*nruns;
Xorg=initialsets[ss,];


y1=0;
do iset=1 to &niset;
as=(iset-1)*&nalts+1:iset*&nalts;
Xlast=Xorg[as,];
call choice(beta,Xlast,y1,iset,y);
y1=y;
end;
cset=&niset;

/*--------------------------------------------------------------*/
/*This following part is only used when there is no initial set.*/
/*(2)*/
/*do j=1 to &niset;
k=(j-1)*&nalts+1:j*&nalts;
ncands = nrowcand/&nalts;
	do h = 1 to ncands until(found);
	l=(h-1)*&nalts+1:h*&nalts;
	found = (Xorg[k,] = holddesign[l,]);
	end;
idx = setdif(1:nrowcand,l);
holddesign = holddesign[idx,];
nrowcand=nrow(holddesign);
end;*/
/*--------------------------------------------------------------*/

do jj=&niset+1 to &tset;

/* The following "logpos", "grd" and "hes" functions are used in Newton-Raphson optimization 
method, called "nlpnrr": see the link for more information:
http://support.sas.com/documentation/cdl/en/imlug/59656/HTML/default/viewer.htm#langref_sect190.htm
The Newton-Raphson optimization is used for estimating the MNL model to update the prior 
repeatedly. */

/*This function is used for computing the loglikelihood of MNL model*/
/*pmean = prior mode*/
start logpos(b) global(Xorg,y,logprior1,pmean,pcov,cset);
logl=0;
	do si = 1 to cset;
    k = (si-1)#&nalts+1 : si#&nalts;  
   	p = exp(Xorg[k,]*t(b));
	p = p/sum(p);
	logl=logl+log(p[y[si,]]);
   	end;
logprior2=-0.5*(b-t(pmean))*inv(pcov)*(t(b)-pmean);
logpost=logl+logprior1+logprior2;
return(logpost);
finish logpos;

/*This function is for computing the gradient - the 1st derivariate of loglikelihood of MNL model
This is required for Newton-Raphson optimization. In case you don't specify the gradient and 
hessian functions in Newton-Raphson, approximations are used, that can cause worse optimization
results. */
start grd(b) global(Xorg,y,pmean,pcov,cset);
nobs=&nalts*cset;
choicematrix=j(nobs,1,0);
do ire=1 to cset;
ses=(ire-1)*&nalts;
ind=ses+y[ire];
 choicematrix[ind,]=1;
end;

   g = j(1, &df, 0);
   utils = Xorg*t(b);
   exputils = exp(utils); 
       	do j = 1 to cset;
	 	som = 0;
     	ssom = 0; 
	 		do k = (j-1)#&nalts+1 to j#&nalts;
	 		som = som + exputils[k];
      		ssom = ssom + exputils[k] * Xorg[k, ];
      		end;
     		do k = (j-1)#&nalts+1 to j#&nalts; 
      		der = choicematrix[k,] * (Xorg[k, ] - ssom/som);
      		g[, ] = g[, ] + der; 
      		end;  
	 	end;
gprior=inv(pcov)*(pmean-t(b));
gpos=(t(g)+gprior)`;
return (gpos); 
finish grd;

/*This function is for computing the hessian - the 2nd derivariate of loglikelihood of MNL model*/
start hes(b) global(Xorg,pcov,cset);

info=j(&df,&df,0);
do scc=1 to cset;
   xx=Xorg[(scc-1)#&nalts+1:scc#&nalts,];
   u=xx*t(b);
   p=exp(u)/sum(exp(u));
   pp=diag(p);
   info=info+(t(xx)*(pp-p*t(p))*xx);
end;
hes=-info-inv(pcov);


return(hes);
finish hes;


/*--------------------------------------------------------------*/
/*In the following part, Importance Sampling procedure is given*/
/*The following parameters are used in Newton-Raphson optimization; they are some options of this
method. */
b=j(1,&df,0); /*initial values are set to zero*/
bmode = j(1, &df, 0.);/*this is the storage matrix for estimates of betas that will be updated prior*/
optn = {1, 0};
tc=repeat(.,1,12);
tc[1]=600;
tc[7]=1.e-8;
call nlpnrr(rc,bmode,"logpos",b,optn,,,,,"grd","hes") tc=tc; /*Compute Beta*_n; Appendix C - (1)in the paper*/ 


 
covst=-inv(hes(bmode)); /*this is -inverse(H), it is the covariance matrix of multivariate
						student t - in other words, covariance of importance density; 
						used in Importance Sampling Appendix C - (1)in the paper*/


prbeta=RANDMVT(&ndraws,&dfv,bmode,covst,z); /*random draws from Multivariate student t with 
											newly estimated betas and -inverse(H) of these 
											estimates; Appendix C - (2)in the paper*/

prbeta=prbeta`;

/*The following 4 lines for creating storage (empty) matrices; required for for-loop results.*/
/*(1)*/denpi0=j(&ndraws,1,0); 
/*(2)*/denlikeli=j(&ndraws,1,0);
/*(3)*/deng=j(&ndraws,1,0);
/*(4)*/wsub=j(&ndraws,1,0);
/*As the the posterior distribution has no close form, the non-parametric way to have the probability
density funcion is used: the density of kernel of q()-the posterior density function*/ 
do ndr=1 to &ndraws; 
denpi0[ndr,]=prior1*exp(-0.5*t(prbeta[,ndr]-pmean)*inv(pcov)*(prbeta[,ndr]-pmean)); /*prior probability; */
denlikeli[ndr,]=prob2(y,Xorg,prbeta[,ndr]); /*likelihood*/
deng[ndr,]=MVST(&dfv,prbeta[,ndr],t(bmode),covst);
wsub[ndr,]=denlikeli[ndr,]*denpi0[ndr,]/deng[ndr,]; /*w_r , Appendix C - (3) in the paper*/
end;
weight=wsub/wsub[+,]; /*w_r , Appendix C - (3) in the paper*/

/*DB-Error computation given in Appendix C - (4) in the paper*/
infomat=j(&ndraws*&df,&df,0);
do dr=1 to &ndraws;
	nrt=(dr-1)*&df+1:dr*&df;
	b=prbeta[,dr];
	call infologit(Xorg,b,cset,info1); /*information matrix of the initial set with priors*/
	infomat[nrt,]=info1;
end;

cset=cset+1;
dsmall=10000;
nsetcand=nrowcand/&nalts;

/*Choose new choice set*/
 do      uu=1 to nsetcand;
		 ij=(uu-1)*&nalts+1:uu*&nalts;
		 X=holddesign[ij,];
	     dmin=0;
	     do dr=1 to &ndraws;
			nrt=(dr-1)*&df+1:dr*&df;
			subinfo=infomat[nrt,];
			b1=prbeta[,dr];
			call infologit2(X,b1,info2);/*information matrix of the (s+1)th set*/
			dmin=dmin+det(info2+subinfo+invcovar)##(-1/&df)*weight[dr,]; /*the DB-error formula*/
		 end;
			if dmin<dsmall then;
		 do;
			dsmall=dmin;
			Xnewsub=X;
			index=ij;
		end;
end;


Xorg=Xorg//Xnewsub; /*the new set is added to initial sets*/
y1=y;
call choice(beta,Xnewsub,y1,cset,y); /*choice data is simulated for the new set*/
idx = setdif(1:nrowcand,index);
holddesign = holddesign[idx,];
nrowcand=nrow(holddesign);
end;


y=y`;
print y;
nxmatsub=&nalts#&tset;
Xmat[(res-1)*nxmatsub+1:res*nxmatsub,]=Xorg;
ymat[res,]=y;
end;

create seqdb.Xseqdes from Xmat; append from Xmat;
create seqdb.Yseqdes from ymat; append from ymat;
options notes;
quit;
