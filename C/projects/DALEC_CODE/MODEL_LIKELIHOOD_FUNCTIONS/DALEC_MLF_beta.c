#pragma once
#include <math.h>
#include "DALEC_ALL_LIKELIHOOD.c"
#include "../../../math_fun/ipow.c"


double DALEC_MLF_beta(DATA DATA,double *PARS){

/*Copy model pointer for brevity*/
DALEC *MODEL=(DALEC *)DATA.MODEL;


/*Setting probability to 0*/
DATA.M_P[0]=0;


/*PARAMETER LOG LIKELIHOOD*/
/*Eventually push to each function(?) OK here for now*/
DATA.M_P[0]=DATA.M_P[0]+LIKELIHOOD_P(DATA,PARS);


/*running model*/
/*For internal checks, i.e. EDCs, are mediated by DATA.M_P[0]*/
DATA.M_P[0]=DATA.M_P[0]+MODEL->dalec(DATA, PARS);

/*Only run likelihood calculation if model is physical*/
if (isinf(DATA.M_P[0])>-1){
/*storing GPP (this should really be done elsewhere)*/
int n;
for (n=0;n<DATA.nodays;n++){DATA.M_GPP[n]=DATA.M_FLUXES[n*DATA.nofluxes];}


/*Calculate likelihood*/
DATA.M_P[0]=DATA.M_P[0]+LIKELIHOOD(DATA);


}

/*Returning the log likelihood P*/
return DATA.M_P[0];


}





/*Not worth spending too much time on, until we get to stacking EDCs*/

double EDC_DALEC_MLF_beta(DATA DATA, double *PARS){

/*Copy model pointer for brevity*/
DALEC *MODEL=(DALEC *)DATA.MODEL;

struct EDCDIAGNOSTIC EDCD;
/*Copy default structure*/
/*EDCD=*((DALEC *)DATA.MODEL)->EDCD;*/
EDCD=*MODEL->EDCD;

/*running model*/
MODEL->dalec(DATA, PARS);

/*LIKELIHOOD (log likelihood)*/
/*EDCs are individually counted*/
/*Only counted if EDCSWITCH is on*/
double P;
int tot_exp=0,n;
for (n=0;n<EDCD.nedc;n++){
tot_exp+=1-ipow(DATA.M_EDCD[n],EDCD.SWITCH[n]);}
P=-0.5*((double)tot_exp*10)*DATA.EDC;

int k;
for (k=0;k<100;k++){printf("%i ",DATA.M_EDCD[k]);}
printf("\n*********\n");


/*overriding if model likelihood is zero or erroneous*/
double ML=DATA.MLF(DATA,PARS);


if (DATA.EDC==0 && ( isinf(ML)==-1 || isinf(ML)==1 || isnan(ML) )){P=P-0.5*10;}
/*if (DATA->EDC==0 && (isinf(ML)==-1 || isnan(ML))){P=P-0.5*10;}
*/

return P;

}




