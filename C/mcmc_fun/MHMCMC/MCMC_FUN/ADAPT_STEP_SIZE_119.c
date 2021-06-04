#pragma once
#include "../../../math_fun/cholesky.c"
#include "../../../math_fun/declare_matrix.c"
int tempfun_print_mat(double **CM, int N){
int p,pp;
for (p=0;p<N;p++){
for (pp=0;pp<N;pp++){
printf("%6.6e ",CM[p][pp]);}
printf("\n");}
printf("\n");
return 0;}

int tempfun_print_vec(double *CM, int N){
int p;
for (p=0;p<N;p++){
printf("%6.6e ",CM[p]);}
printf("\n");
return 0;}


int statfun_incremental_cov_mat(double **CM, int Npars, int Nsamples, double *CMmean, double *newsample, double ar){


/*Step 1. Calculate new mean*/
int p,pp;
double N=(double)Nsamples;
/*Mi=(M*(N)+ ar*x)/(N+ar);*/
double *CMmeannew=calloc(Npars, sizeof(double));

for (p=0;p<Npars;p++){
CMmeannew[p]=(CMmean[p]*N + ar*newsample[p])/(N + ar);}

/*Step 2. calculate new cov mat*/
/*CMOUT=CM*(N-1)/(N-1+ar)+((N)* M'*M- (N+ar)*Mi'*Mi + x'*x*ar)/(N-1+ar);*/
for (p=0;p<Npars;p++){
for (pp=0;pp<Npars;pp++){
CM[p][pp]=CM[p][pp]*(N-1)/(N-1+ar) +(N*CMmean[p]*CMmean[pp] - (N+ar)*CMmeannew[p]*CMmeannew[pp] + ar*newsample[p]*newsample[pp])/(N-1+ar);}}


/*Update CMmean and clear CMmeannew*/
for (p=0;p<Npars;p++){CMmean[p]=CMmeannew[p];}
free(CMmeannew);

return 0;}



int ADAPT_STEP_SIZE_119(double *PARSALL, PARAMETER_INFO PI, COUNTERS * N, MCMC_OPTIONS MCO){

/*first adjusting all step sizes*/

int n,p,pp;
double fac=1, scfac=(1-(double)N->ACC/(double)MCO.nOUT/2);
double adapfac=scfac*0.001+1;
/*ratio loc vs total*/
/*
double adapfac=1.5;*/


/*MAKE THIS COMPATIBLE WITH STRUCTS AND CODE*/



/*Scale factor adaptation*/
if (N->ACCLOC>PI.npars*10){
	/*Growth in all dims if ACCRATE>0.23*/
	if (N->ACCRATE<0.23){N->amp=N->amp/adapfac;}
	else {N->amp=N->amp*adapfac;}}
/*Done with scale factor adaptation*/
N->amp=1;



/*Only derive covariance once*/
if (N->ACCLOC>PI.npars*10 & N->covsamplecount==0){

	double **norparvec;
	declare_matrix(&norparvec,PI.npars, N->ACCLOC);
	int i, os;
	/*Offset (gaining 2, losing 1)*/
	os=(int)floor((double)N->ACCLOC/2);
	/*std of each parameter (for all recently N->ACCLOCepted parameters)*/
	for (p=0;p<PI.npars;p++){
	for (n=0;n<(N->ACCLOC-os);n++){
	/*transforming parameters and storing in pointer*/
	norparvec[p][n]=par2nor(PARSALL[PI.npars*(n+os)+p],PI.parmin[p],PI.parmax[p]);}}

	/*Calculare covariance of norparvec AND store*/
	for (p=0;p<PI.npars;p++){
	for (pp=0;pp<PI.npars;pp++){
	N->parcovariance[p][pp]=covariance(norparvec[p],norparvec[pp],N->ACCLOC-os);}}


	/*Calculate and store mean*/
	for (p=0;p<PI.npars;p++){
        for (n=0;n<(N->ACCLOC-os);n++){
	N->parmean[p]=mean(norparvec[p],N->ACCLOC-os);}}
	N->covsamplecount=N->ACCLOC-os;

	free(*(norparvec));
	free(norparvec);
}
/*Done with Covariance derivation*/



/*Now only deriving incremental covariance changes*/
/*statfun_incremental_cov_mat(double **CM, int Npars, int Nsamples, double *CMmean, double *newsample, double ar)*/

else if (N->ACCLOC>PI.npars*10 & N->covsamplecount>0){
	


	/*Derive terms needed for statfun*/
	/*Normalized parameter vector sample to add/remove*/
	double *norparveci=calloc(PI.npars,sizeof(double));
	/*add and (maybe) remove sample*/
	int isample=N->ACCLOC;
	for (p=0;p<PI.npars;p++){norparveci[p]=par2nor(PARSALL[PI.npars*(isample-1)+p],PI.parmin[p],PI.parmax[p]);}

	
	/*Step 3. Use new sample to update covariance matrix*/
	/*Note that parmean and parcovariance will get updated here*/
        
	/*SAVE A COPY OF CM*/


	statfun_incremental_cov_mat(N->parcovariance,PI.npars,N->covsamplecount,N->parmean,norparveci,1);
	N->covsamplecount=N->covsamplecount+1;

	/*Step 4. Remove parameter if mod(N,2) is not zero*/
	if (N->ACCLOC % 2>0){
	/*Copying matrix*/	

	isample=(int)floor((double)N->ACCLOC/2);
	for (p=0;p<PI.npars;p++){norparveci[p]=par2nor(PARSALL[PI.npars*(isample-1)+p],PI.parmin[p],PI.parmax[p]);}
	statfun_incremental_cov_mat(N->parcovariance,PI.npars,N->covsamplecount,N->parmean,norparveci,-1);
        N->covsamplecount=N->covsamplecount-1;
	



	}


}

/*
if (N->ACCLOC>10*PI.npars){
printf("N->ACCLOC = %i\n",N->ACCLOC);
printf(" N->parcovariance = %3.3f\n",N->parcovariance[0][0]);
}*/
/*Now calculating cholesky and storing*/

if (N->ACCLOC>PI.npars*10){

double **CCM;
declare_matrix(&CCM,PI.npars,PI.npars);
int choldiag=cholesky(N->parcovariance,PI.npars,CCM);

if (choldiag==1){
for (p=0;p<PI.npars;p++){
for (pp=0;pp<PI.npars;pp++){
N->parcholesky[p][pp]=CCM[p][pp];}}}


	free(*(CCM));free(CCM);
}
/*Done with Cholesky calculation*/


return 0;

}





