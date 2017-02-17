libname klchoice 'X:\KBI_ORSTAT\Marjolein\KL for Choice\Programs';

%let nalts = 2;
%let nsets = 5;

%let nres = 1;

%let df = 3;

%let nrandes = 1; 
%let ndraws = 512;


proc plan ordered;
factors x1=5 x2=4 x3=4/noprint;
output out=candidate;
run;

proc print data= candidata; 
run;

data candidate; set candidate;
if x1=1 then x1=10; if x1=2 then x1=50; if x1=3 then x1=90; if x1=4 then x1=130; if x1=5 then x1=170; 
if x2=1 then x2=10; if x2=2 then x2=40; if x2=3 then x2=70; if x2=4 then x2=100;
if x3=1 then x3=log(0.5); if x3=2 then x3=log(1.5); if x3=3 then x3=log(2.5); if x3=4 then x3=log(3.5);
run;

proc iml;

*use candid(keep = &_trgind); 
use candidate;
read all into candidat;
use klchoice.latfull_3par; read all into zfull;


use klchoice.mutrue3_new; read all into pmean;
use klchoice.sigmatrue3_new; read all into pcov;


nobs = &nsets#&nalts;

do i = 1 to nrow(candidat)-1;
	do j = 1 to nrow(candidat)-i;
    	holddesignfull = holddesignfull // candidat[i, ];
    	holddesignfull = holddesignfull // candidat[i + j, ];
	end;
end;

ncandshold = nrow(holddesignfull)/&nalts;

start Kl(kl,design,beta);
nset=nrow(design)/&nalts;
kl = 0;  
do s = 1 to nset;
	k = (s-1)#&nalts+1 : s#&nalts;
	X = design[k,];
	numer = X*beta;
	mmat = numer[<>,];
	numermax = exp(numer - mmat);
	denom = numermax[+,];
	probs = numermax/denom;
	logprobs = log(probs);
	elemfin = probs[,:];
	logelemfin = logprobs[,:];
	klsub = 0;
	do klelem = 1 to &nalts;
	klsub = klsub + ( elemfin[klelem,] * (log(elemfin[klelem,]) - logelemfin[klelem,]) );
	end;
kl = kl + klsub;
end;
finish Kl;

start finddesign(des,indvec,holddesignfull);
r = (indvec-1)*&nalts+1 || indvec*&nalts;
*r = (indvec-1)*&nalts+1 || (indvec-1)*&nalts+2 || indvec*&nalts;
indices = shape(r,&nsets*&nalts,1);
des = holddesignfull[indices,];
finish finddesign;

design_total = j(&nres*&nsets*&nalts,&df,.);

do res = 1 to &nres;
klmin = -10000;
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
	call Kl(kl,randes,beta);
  	oldkl = kl; 
   		do iteration = 1 to 7 until(converge);	  
    		do desi = 1 to &nsets;
	 		besttry = indvec[desi, ]; 
     			do candidi = 1 to ncandshold; 
					find = (indvec=candidi);
					ind = find[+,];
					if ind=0 then do;
	  					indvec[desi, ] = candidi;
						call finddesign(subdes,indvec,holddesignfull);
						call Kl(trykl,subdes,beta);	
	  					if trykl > kl then do; 
	   					kl = trykl;
	   					besttry = candidi;
	   					end; 
      				end; 
				end; 
	 		indvec[desi, ] = besttry;
	 		end;
		currentkl = kl; 
    	converge =((currentkl - oldkl) / max(currentkl, 1e-8) < 0.005);
		oldkl = currentkl;
		end;
    	if kl > klmin then do;
		call finddesign(des,indvec,holddesignfull);
    	design_total[((res-1)*&nsets*&nalts)+1:res*&nsets*&nalts,] = des;
    	klmin = kl;
		end;
	end;
end;
print design_total;

create klchoice.initialsets_KLP from design_total;
print klchoice.initialsets_KLP;
append from design_total;

quit;

