#pragma once
#include <math.h>
#include "../../../math_fun/randn.c"

int STEP_119(double *pars0, double *pars, PARAMETER_INFO PI,COUNTERS N,  MCMC_OPTIONS MCO){

/*FIXEDPARS*/
/*ones and zeros depending on whether parameters are kept fixed*/

int n,nn,fp,withinlim=1;
double npar,*rn,*rn0,step, npar_new, rn2;


rn0=calloc(PI.npars,sizeof(double));
rn=calloc(PI.npars,sizeof(double));


/*cholesky sampling here*/
/*cholesky matrix = upper triangular*/
/*correlated random numbers = rn * U*/

for (n=0;n<PI.npars;n++){rn0[n]=randn();}

/*looping through each rn parameter*/
for (n=0;n<PI.npars;n++){
/*multiply each row by chol column*/
for (nn=0;nn<n+1;nn++){
rn[n]+=rn0[nn]*N.parcholesky[nn][n];}}


/*
if (N.parcholesky[0][0]==12345){

for (nn=0;nn<PI.npars;nn++){
for (n=0;n<PI.npars;n++){
printf("%3.3f ",N.parcholesky[nn][n]*1e4);}
printf("\n");}
printf("***chol**\n");
for (nn=0;nn<PI.npars;nn++){
for (n=0;n<PI.npars;n++){
printf("%3.3f ",N.parcovariance[nn][n]*1e4);}
printf("\n");}
printf("**cov***\n");
for (n=0;n<PI.npars;n++){printf("%3.3f ",rn[n]*1e4);};printf("\n**rn***\n");
for (n=0;n<PI.npars;n++){printf("%3.3f ",rn0[n]*1e4);};printf("\n**rn0***\n");
}
*/

double sd=pow(2.38,2)/(double)PI.npars;

/*SAMPLING PARAMETERS*/
for (n=0;n<PI.npars;n++){
	/*sampling for each parameter*/
	/*normalizing parameter*/
		rn2=randn();
		step = rn[n]*N.amp*sd + rn2*MCO.minstepsize;
		npar=par2nor(pars0[n],PI.parmin[n],PI.parmax[n]);
		/*circular parameters (e.g. day of year)*/
		npar_new=step+npar;
		if (npar_new<0 | npar_new>1){withinlim=0;}
/*printf("rn/rn2/npar/step/new = %2.2f/%2.2f/%2.2f/%2.2f/%2.2f\n",rn[n],rn2,npar,step,npar_new);*/
		pars[n]=nor2par(npar_new,PI.parmin[n],PI.parmax[n]);}

free(rn);free(rn0);
return withinlim;


}


