#pragma once
#include "../../../DALEC_CODE/DALEC_ALL/mean_pool.c"
#include "../../../DALEC_CODE/DALEC_ALL/mean_annual_pool.c"
#include "../../../DALEC_CODE/DALEC_ALL/expdecay2.c"
#include "../../../../math_fun/std.c"
#include "../../../../math_fun/ipow.c"
#include "stdlib.h"
#include "stdio.h"


int EDC2_1005(double const *pars, DATA DATA, struct EDCDIAGNOSTIC *EDCD)
{

/*Extract DALEC model here*/
/*Copy model pointer for brevity*/
DALEC *MODEL=(DALEC *)DATA.MODEL;

double *MET=DATA.MET;
double *POOLS=DATA.M_POOLS;
double *FLUXES=DATA.M_FLUXES;
int nodays=DATA.nodays;
double *parmax=DATA.parmax;
double meantemp=DATA.meantemp;

/*EDCD=EDCD2;*/

/* EDC2_SUMMARY
6. Cfol:Croo
7-14. Pool trajectory EDCs
15. Ensuring that wilting point is at or below the mean H2O pool
35.
*/

/*THESE EDC2 checks are for DALEC_FIREBUCKET3d*/
int EDC=1,n=0,m=0,edc=0;
int DIAG=EDCD->DIAG;/*1 or 0*/



/*FIREBUCKET*/
int nomet=MODEL->nomet;
int nopools=MODEL->nopools;
int nofluxes=MODEL->nofluxes;
int done=0;
int k=0;



/*deriving mean pools here!*/
double *MPOOLS;
MPOOLS=calloc(nopools,sizeof(double));
if (MPOOLS==0){printf("WARNING NULL POINTER");}
for (n=0;n<nopools;n++){MPOOLS[n]=mean_pool(POOLS,n,nodays+1,nopools);};

/*deriving mean January pools*/
/*Assuming COMPLETE years*/
double *MPOOLSjan;
/*pool interval*/
int dint=(int)floor(nodays/(MET[nomet*(nodays-1)]-MET[0])*365.25);
/*declaring mean pool array*/
MPOOLSjan=calloc(nopools,sizeof(double));if (MPOOLSjan==0){printf("WARNING NULL POINTER");}
/*deriving mean jan pools*/
/*based on all jan pools except initial conditions*/
for (n=0;n<nopools;n++){
for (m=0;m<(nodays/dint+1);m++){
MPOOLSjan[n]=MPOOLSjan[n]+POOLS[nopools*(m*dint)+n]/(nodays/dint+1);}}
/*printing just to make sure*/
/*for (n=0;n<nopools;n++){printf("Pool = %d, janmnean=%2.2f\n",n,MPOOLSjan[n]);}*/


/***********************EDCs start here****************************/


/*EDC no 6*/
/*0.2*Cf < Cr < 5*Cf*/
/*Cfoliar : Croot = 5:1 or 1:5*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[6-1]==1)) & (MPOOLS[1]>MPOOLS[2]*5 | MPOOLS[1]*5<MPOOLS[2])){EDC=ipow(0,EDCD->SWITCH[6-1]);EDCD->PASSFAIL[6-1]=0; }


/*equilibrium factor (in comparison to C_initial)*/
double EQF=EDCD->EQF;

/*Total fluxes*/
double *FT;
FT=calloc(nofluxes,sizeof(double));
int f=0;
for (f=0;f<nofluxes;f++){FT[f]=0;for (n=0;n<nodays;n++){FT[f]+=FLUXES[n*nofluxes+f];}}
/*Total prec*/
double PREC=0;
for (n=0;n<nodays;n++){PREC+=MET[n*nomet+8];}


double Fin[8];
double Fout[8];
double Pstart;
double Pend;
/*temporary print switch*/
int psw=0;
/*exponential decay tolerance*/
double etol=0.1;





/*Inputs and outputs for each pool*/
/*labile*/
Fin[0]=FT[4];
Fout[0]=FT[7]+FT[17]+FT[23];
/*foliar*/
Fin[1]=FT[3]+FT[7];
Fout[1]=FT[9]+FT[18]+FT[24];
/*root*/
Fin[2]=FT[5];
Fout[2]=FT[11]+FT[19]+FT[25];
/*wood*/
Fin[3]=FT[6];
Fout[3]=FT[10]+FT[20]+FT[26];
/*litter*/
Fin[4]=FT[9]+FT[11]+FT[23]+FT[24]+FT[25];
Fout[4]=FT[12]+FT[14]+FT[21]+FT[27];
/*som*/
Fin[5]=FT[10]+FT[14]+FT[26]+FT[27];
Fout[5]=FT[13]+FT[22];
/*PAH2O*/
Fin[6]=PREC;
Fout[6]=FT[28]+FT[29]+FT[30];
/*PUH2O*/
Fin[7]=FT[30];
Fout[7]=FT[31];


/*Inlcuding H2O pool*/
/*EDCs 7-13 - inputs, outputs and exponential tolerance*/

/*mean input/output ratio and start ratio*/
double Rm, Rs;



for (n=0;n<nopools;n++){
/*start and end pools*/
Pstart=POOLS[n];
Pend=POOLS[nopools*nodays+n];
/*mean input/output*/
Rm=Fin[n]/Fout[n];
/*Theoretical starting input/output*/
Rs=Rm*MPOOLSjan[n]/Pstart;

if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[7-1+n]==1))
& ((fabs(log(Rs))>log(EQF)) || (fabs(Rs-Rm)>etol)))
{EDC=ipow(0,EDCD->SWITCH[7-1+n]);EDCD->PASSFAIL[7-1+n]=0;}


/*storing EDCPROB: i.e. the log probability of each EDC based on a gaussian representation*/
/*of each constraint*/
EDCD->EDCPROB[7-1+n]=-0.5*pow(log(Rs)/log(EQF),2);/*-0.5*pow((Rs-Rm)/etol,2);*/

if (psw==1){
printf("****\n");
printf("Pool %i EDCDPROB = %f\n",n,EDCD->EDCPROB[7-1+n]);
printf("Pool %i Fin = %f,Fout = %f\n",n+1,Fin[n],Fout[n]);
printf("Pool %i Pstart = %f,Pend = %f\n, Pmeanjan=%f\n",n+1,Pstart,Pend,MPOOLSjan[n]);
printf("fabs(log(Fin/Fout)) = %f\n",fabs(log(Fin[n]/Fout[n])));
printf("fabs(log(Pend/Pstart)) = %f\n",fabs(log(Pend/Pstart)));
printf("Rm = %f\n",Rm);
printf("Rs = %f\n",Rs);
printf("log(EQF) = %f\n",log(EQF));
printf("etol = %f\n",etol);
printf("****\n");}}


/*Ensuring that wilting point is at or below the mean H2O pool EDC14*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[15-1]==1)) & (pars[25]>MPOOLS[6])){EDC=ipow(0,EDCD->SWITCH[15-1]);EDCD->PASSFAIL[15-1]=0;}

/***********************EDCs done here****************************/



/*Additional faults can be stored in positions 35-40*/

/*PRIOR RANGES - ALL POOLS MUST CONFORM*/
int pidx[]={17,18,19,20,21,22,26,35};

for (n=0;n<nopools-1;n++){if ((EDC==1 || DIAG==1) & ((MPOOLS[n])>parmax[pidx[n]])){EDC=0;EDCD->PASSFAIL[35-1]=0;}}


int PEDC;
/*ensuring minimum of each pool is zero & finite*/
if (EDC==1 || DIAG==1)
{double min; int nn;n=0;
while ((n<nopools) & (EDC==1 || DIAG==1))
{nn=0;PEDC=1;while ((nn<nodays+1) & (PEDC==1))
{if ((POOLS[n+nn*nopools]<0) || isnan(POOLS[n+nn*nopools])==1)
{EDC=0;PEDC=0;EDCD->PASSFAIL[35+n]=0;}nn=nn+1;};
n=n+1;
}
}
/*
if (EDC==1){printf("\n");oksofar("EDC(36-43) passed");}
*/








/*FREE MEMORY*/
free(FT);
free(MPOOLS);
free(MPOOLSjan);





/*final check confirming EDC = 1 or 0*/
int Num_EDC=100;
if (DIAG==1){for (n=0;n<Num_EDC;n++){if (EDCD->PASSFAIL[n]==0){EDC=0;}}}




/*Returning EDC */
return EDC;

}







