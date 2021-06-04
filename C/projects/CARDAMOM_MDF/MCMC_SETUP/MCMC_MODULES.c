

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CARDAMOM_GENERAL/CARDAMOM_READ_BINARY_DATA.c"
#include "PROJECT_FUN/FIND_EDC_INITIAL_VALUES.c"
#include "../../../math_fun/declare_matrix.c"
/*DALEC_SYNTHETIC SETUP*/


int READ_PARI_DATA(PARAMETER_INFO *PI, DATA *CARDADATA,MCMC_OUTPUT *MCOUT,MCMC_OPTIONS *MCOPT, char *CLA[]){
/*READING IN DALEC_SYNTHETIC CARDADATA*/
/*opening file*/


/*CHANGE4EDDIE*/
/*this only applies to the native environment - change accordingly to add default file!*/
char filename[200];
if (atoi(CLA[0])<1){strcpy(filename,"MCMC_SETUP/TEST_BINARY_DATASET.cbf");}
else{strcpy(filename,CLA[1]);}





/*defining initial values*
 * Need to perform MCMC run to determine this*/
FIND_EDC_INITIAL_VALUES(*CARDADATA,PI,MCOPT);


/*resetting PI-stepsize (as this has been changed)*/
int n;
for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.0001;}

oksofar("Done with initial parameters");


/*This function is in the MCMC folder*/
INITIALIZE_MCMC_OUTPUT(*PI,MCOUT,*MCOPT);

for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.0001;}

return 0;
}


/*MCMC OPTIONS*/
int READ_MCOPT(MCMC_OPTIONS *MCOPT, char *CLA[]){
/*number of command line imputs*/
int ncli=atoi(CLA[0]);

/*defining MCMC_OPTIONS structure*/
MCOPT->APPEND=0;
MCOPT->nADAPT=100;
MCOPT->fADAPT=0.5;
/*command line (or default) values*/
if (ncli<3){MCOPT->nOUT=1000;}else{MCOPT->nOUT=atoi(CLA[3]);};
if (ncli<4){MCOPT->nPRINT=1000;}else{MCOPT->nPRINT=atoi(CLA[4]);};
if (ncli<5){MCOPT->nWRITE=10;}else{MCOPT->nWRITE=atoi(CLA[5]);};
if (ncli<6){MCOPT->minstepsize=1e-7;}else{MCOPT->minstepsize=atof(CLA[6]);};
if (ncli<7){MCOPT->mcmcid=1;}else{MCOPT->mcmcid=atof(CLA[7]);};
if (ncli<8){MCOPT->nADAPT=100;}else{MCOPT->nADAPT=atoi(CLA[8]);};

MCOPT->randparini=0;
MCOPT->returnpars=0;
MCOPT->fixedpars=0;
char outfile[200], stepfile[200],startfile[200];
if (ncli<2){strcpy(outfile,"MOUT_");strcpy(stepfile,"MOUT_");}
else{strcpy(outfile,CLA[2]); strcpy(stepfile,CLA[2]); strcpy(startfile,CLA[2]);}
/*STEP file comes with SUFFIX*/
/*PARS file no longer does!*/
/*strcat(outfile,"PARS");*/
strcat(stepfile,"STEP");
strcat(startfile,"START");
/*directory*/
strcpy(MCOPT->outfile,outfile);
strcpy(MCOPT->stepfile,stepfile);
strcpy(MCOPT->startfile,startfile);

return 0;

}


/*Enter all fields originally defined with MALLOC*/
int MEMORY_CLEANUP(DATA CARDADATA, PARAMETER_INFO PI, MCMC_OPTIONS MCOPT, MCMC_OUTPUT MCOUT){

free(PI.parmin);
free(PI.parmax);
free(PI.parini);
free(PI.parfix);
free(PI.stepsize);
free(PI.transform);

FREE_DATA_STRUCT(CARDADATA);
/*
free(CARDADATA.MET);
free(CARDADATA.LAI);
free(CARDADATA.NEE);
free(CARDADATA.WOO);
free(CARDADATA.GPP);


free(CARDADATA.M_FLUXES);
free(CARDADATA.M_LAI);
free(CARDADATA.M_NEE);
free(CARDADATA.M_POOLS);
free(CARDADATA.M_GPP);

free(CARDADATA.parmin);
free(CARDADATA.parmax);

if (CARDADATA.ngpp>0){free(CARDADATA.gpppts);}
if (CARDADATA.nlai>0){free(CARDADATA.laipts);}
if (CARDADATA.nnee>0){free(CARDADATA.neepts);}
if (CARDADATA.nwoo>0){free(CARDADATA.woopts);}

*/
free(MCOUT.best_pars);


return 0;}




/*this function initializes the PI fields
It is called from DALEC_ALL_TEMPLATE or equivalent higher level function*/
int INITIALIZE_PI_STRUCT(PARAMETER_INFO * PI, DATA * CARDADATA, MCMC_OPTIONS *MCO){
oksofar("initializing PI stucture");
/*contains 6 fields with min max log for par and par*/
PI->parmin=calloc(CARDADATA->nopars,sizeof(double));
PI->parmax=calloc(CARDADATA->nopars,sizeof(double));
PI->parini=calloc(CARDADATA->nopars*MCO->nchains,sizeof(double));
PI->parfix=calloc(CARDADATA->nopars,sizeof(double));
PI->stepsize=calloc(CARDADATA->nopars,sizeof(double));
PI->transform=calloc(CARDADATA->nopars,sizeof(int));
/*MAtrix double-pointer allocation*/
oksofar("about to declare matrix");
oksofar("Just declared matrix");


/*defining step size*/
oksofar("fields declared");

int n,m;
PI->npars=CARDADATA->nopars;
for (n=0;n<CARDADATA->nopars;n++){
/*copying minimum and maximum parameter values from CARDADATA*/
PI->parmin[n]=CARDADATA->parmin[n];
PI->parmax[n]=CARDADATA->parmax[n];
PI->stepsize[n]=0.01;}

/*Try first hardcoding transform*/
/*PI->transform[11]=2;
PI->transform[14]=2;
*/


oksofar("stepsizes initialized");
return 0;}



