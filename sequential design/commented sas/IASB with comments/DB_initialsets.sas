libname seqdb 'C:\Users\u0105757\Desktop\sequential design\commented sas\IASB with comments';

%let nalts = 2; /* # of alternatives*/
%let nsets = 5; /* # of initial choice sets*/

%let nres = 1; /* # of respondents; of course 10 respondents won't be enough tohave accurate 
estimates; mostly it is chosen more than 50; as 200 or 250*/ 

%let df = 4; /* # of betas*/

%let nrandes = 1; /* # of random starts (candidate designs) used for the Modified Federov Algorithm;
It is mostly set to 1000 or for quick studies, it is stated as 300 or 100 */
%let ndraws = 16; /* # of random parameter draws. For lattice points, this number is always 
set to 2^m, here m=4, mostly it is used as m=8 or 9 */

/*This following part generates the full-factorial design for the specified design setting
For a small design set-up, 3^2/2/10 is used. 2 attributes with 2 levels, 2 alternatives and 
5 choice sets.*/
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
proc iml;
/* The effect coded full-factoral design is used for creating a whole (possible combination of
choice sets)candidate sets  */
use candid(keep = &_trgind); read all into candidat;
/*zfull = all lattice points for random draws*/
use seqdb.latfull; read all into zfull;

use seqdb.mutrue; read all into pmean; /*prior means used for generating design*/
use seqdb.sigmatrue; read all into pcov; /*prior covariance matrix used for generating design*/

nobs = &nsets#&nalts;
/*holddesignfull is the whole (possible combination of choice sets)candidate sets  */ 
do i = 1 to nrow(candidat)-1;
	do j = 1 to nrow(candidat)-i;
    	holddesignfull = holddesignfull // candidat[i, ];
    	holddesignfull = holddesignfull // candidat[i + j, ];
	end;
end;

ncandshold = nrow(holddesignfull)/&nalts;

/*In SAS/IML, the functions starts with "start function_name(inputs output)", ends with "finish"*/
/*This function is used for calculating the D-error of MNL/CL (multinomial logit or conditional
logit) model*/
start Der(derror,design,beta,pcov);
nset=nrow(design)/&nalts;
sumderror = 0;  
np = 1/ncol(design);
do d = 1 to &ndraws;
infomat = j(ncol(design), ncol(design), 0);
	do s = 1 to nset;
	k = (s-1)#&nalts+1 : s#&nalts;
	p = exp(design[k,]*beta[, d]);
	p = p/sum(p);
	z = design[k, ];
	infoset = z`*(diag(p)-p*p`)*z;
	infomat = infomat+infoset; 
	end;
	infomat = infomat + inv(pcov);
    detinfomat = det(infomat);
    dbeff = detinfomat ## np;
    if dbeff <= 0 then dberror = .; else dberror = 1 / dbeff;
    sumderror = sumderror + dberror;    
end;
derror = sumderror / &ndraws;
finish Der;

/* This function is used for combining randomly selected choice sets from holddesignfull
and forming a design */
start finddesign(des,indvec,holddesignfull);
r = (indvec-1)*&nalts+1 || indvec*&nalts;
*r = (indvec-1)*&nalts+1 || (indvec-1)*&nalts+2 || indvec*&nalts;
indices = shape(r,&nsets*&nalts,1);
des = holddesignfull[indices,];
finish finddesign;

design_total = j(&nres*&nsets*&nalts,&df,.);

/*MODIFIED FEDEROV ALGORITHM*/
do res = 1 to &nres;
dmin = 10000;
k = (res-1)*&ndraws+1:res*&ndraws;
z = zfull[k,];
z=z`;
chol = t(root(pcov));
beta = pmean + (chol*z);
	do desnum = 1 to &nrandes;
  	do until (ind=0);
	indvec = ceil(ncandshold*uniform(j(&nsets,1,0)));
	ind=0;
	do i=1 to (&nsets-1);
		do j=i+1 to &nsets;
			if indvec[i,]=indvec[j,] then ind=1;
		end;
	end;
	end;
	call finddesign(randes,indvec,holddesignfull);
	call Der(err,randes,beta,pcov);
  	olderror = err; 
   		do iteration = 1 to 7 until(converge);	  
    		do desi = 1 to &nsets;
	 		besttry = indvec[desi, ]; 
     			do candidi = 1 to ncandshold; 
					find = (indvec=candidi);
					ind = find[+,];
					if ind=0 then do;
	  					indvec[desi, ] = candidi;
						call finddesign(subdes,indvec,holddesignfull);
						call Der(tryerr,subdes,beta,pcov);	
	  					if tryerr < err then do; 
	   					err = tryerr;
	   					besttry = candidi;
	   					end; 
					end;
      			end;     
	 		indvec[desi, ] = besttry;
	 		end;
		currenterror = err; 
    	converge =((olderror - currenterror) / max(olderror, 1e-8) < 0.005);
		olderror = currenterror;
		end;
    	if err < dmin then do;
		call finddesign(des,indvec,holddesignfull);
    	design_total[((res-1)*&nsets*&nalts)+1:res*&nsets*&nalts,] = des;
    	dmin = err;
		end;
	end;
end;

/*This writes out the generated intial choice sets*/
create seqdb.initialsets_DB from design_total;
append from design_total;

quit;

/*Note: In SAS, for-loops occur as do-end loops. Every "do" loop - ends with "end"*/
/*Inside the "proc iml", "print" can be used to check the outputs*/
