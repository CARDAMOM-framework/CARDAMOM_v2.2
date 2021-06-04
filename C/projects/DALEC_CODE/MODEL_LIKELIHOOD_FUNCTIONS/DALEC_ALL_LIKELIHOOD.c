#pragma once
#include "../../../math_fun/std.c"
#include "../../../math_fun/mean.c"
#include "../../../math_fun/max.c"
#include "DALEC_806_MFCF.c"
#include "DALEC_807_MFCF.c"
#include "DALEC_LIKELIHOOD_GRACE_EWT.c"
#include "DALEC_LIKELIHOOD_ET.c"
#include "DALEC_LIKELIHOOD_GPP.c"
#include "DALEC_LIKELIHOOD_LAI.c"
#include "DALEC_LIKELIHOOD_NEE.c"
#include "DALEC_LIKELIHOOD_CH4.c"
/*Any likelihood functions used in multiple MLF functions are kept here!*/



int EDCSETUP(DATA DATA, struct EDCDIAGNOSTIC ** EDCD){

/*TO DO HERE: set EDCD.EDCSWITCHES to DATA.OTHERPRIORS fields 30 to 50*/
/*DON'T FORGET TO DO THE SAME IN CDEA files*/
/*see if it is possible to initiate both using the same function*/
/*E.g. INITIATE_EDC_SWITCHES*/

static struct EDCDIAGNOSTIC EDCDmem;
EDCDmem.DIAG=DATA.EDCDIAG;
int n;
/*EDC switches are stored in DATA->parpriors (positions 31-50);*/
/*switch all on if EDCDIAG<2 , otherwise, switch on/off according to EDCSWITCH*/
/*EDCSETUP function will be included in DALEC_ALL_LIKELIHOOD.c*/
for (n=0;n<100;n++){EDCDmem.SWITCH[n]=1;}
if (EDCDmem.DIAG==2){for (n=0;n<20;n++){EDCDmem.SWITCH[n]=DATA.otherpriors[n+30];}}

/*EQF is stored in the "DATA.otherpriorunc" fields associated with EDCs 7-12*/
/*default value is 2*/
EDCDmem.EQF=DATA.otherpriorunc[36]; if (EDCDmem.EQF==-9999){EDCDmem.EQF=2;}

printf("EDCD->EQF = %2.2f\n",EDCDmem.EQF);

/*Default structure has all EDC=1, and 100 EDCs*/
EDCDmem.nedc=100;
for (n=0;n<EDCDmem.nedc;n++){EDCDmem.PASSFAIL[n]=1;}


/*Double pointer = contents of EDCD reassigned to static struct*/
*EDCD=&EDCDmem;


return 0;
}


double LIKELIHOOD_P(DATA DATA,double *PARS)
{
/*remember - LOG likelihood*/
double p=0,p_lma,pn;
int n;

/*looping through all priors for P*/
/*where no prior distribution is used, insert 9999*/
for (n=0;n<50;n++){if (DATA.parpriors[n]>-9999 & DATA.parpriorunc[n]!=-9999){
if (DATA.parpriorunc[n]>0){
/*log-normal if uncertainty value is positive*/
p=p-0.5*pow(log(PARS[n]/DATA.parpriors[n])/log(DATA.parpriorunc[n]),2);}
else {
/*log-normal if uncertainty value is positive*/
p=p-0.5*pow((PARS[n]-DATA.parpriors[n])/(DATA.parpriorunc[n]),2);}

}}

/*for any other priors, explicitly define functions based on values in DATA.otherpriors*/

/*total biomass (pools 1:4) defined here - using first space of pappriors for total biomass*/
if (DATA.otherpriors[0]>-9999){
p=p-0.5*pow(log((PARS[17]+PARS[18]+PARS[19]+PARS[20])/DATA.otherpriors[0])/log(DATA.otherpriorunc[0]),2);}

return p;}


double LIKELIHOOD(DATA D)
{
int n,dn,m,f;
double P=0, Ps;
double tot_exp;
double CPG;
double mfire=0;
double msif=0;
double mgpp=0;
double mnpdf=0;
double mlai=0;
int npdfi[7]={9,10,11,23,24,25,26};



P=P+DALEC_LIKELIHOOD_GPP(D);
P=P+DALEC_LIKELIHOOD_LAI(D);
P=P+DALEC_LIKELIHOOD_ET(D);
P=P+DALEC_LIKELIHOOD_NEE(D);
if (D.ID==1010){
	P=P+DALEC_LIKELIHOOD_CH4(D);
	}

	if (D.ID==1011){
	P=P+DALEC_LIKELIHOOD_CH4(D);
	}
/*shuang: 101010 was created for climate sensitivity test Nov2020*/

double mam=0,am=0;

/*printf("ETprob = %2.2f\n",P);*/


/*printf("NEEprob = %2.2f\n",P);*/

/*Cbiomass likelyhood*/
/*(previously just woody pool)*/
double biomass;
tot_exp=0;
if (D.nwoo>0){
/*Looping through all available constraints*/
for (n=0;n<D.nwoo;n++){dn=D.woopts[n];
/*biomass = totals of labile, foliar, root, wood*/
biomass=D.M_POOLS[D.nopools*(dn+1)+0]+D.M_POOLS[D.nopools*(dn+1)+1]+D.M_POOLS[D.nopools*(dn+1)+2]+D.M_POOLS[D.nopools*(dn+1)+3];
/*Model-data mismatch*/
tot_exp+=pow(log(biomass/D.WOO[dn])/log(D.otherpriorunc[1]),2);}
//tot_exp+=pow(log(max(biomass,0.01)/max(D.WOO[dn],0.01))/log(max(D.otherpriorunc[1],0.1)),2);} /*shuang for correcting -inf*/
/*adding cost to overall probability estimate*/
P=P-0.5*tot_exp;
}


double som;
tot_exp=0;
if (D.nsom>0){
/*Looping through all available constraints*/
for (n=0;n<D.nsom;n++){dn=D.sompts[n];
/*biomass = totals of litter and som*/
som=D.M_POOLS[D.nopools*(dn+1)+4]+D.M_POOLS[D.nopools*(dn+1)+5];
/*Model-data mismatch*/
tot_exp+=pow(log(som/D.SOM[dn])/log(D.otherpriorunc[7]),2);}
//tot_exp+=pow(log(max(som,0.1)/max(D.SOM[dn],0.1))/log(max(D.otherpriorunc[7],100)),2);} /*shuang for correcting -inf*/
/*adding cost to overall probability estimate*/
P=P-0.5*tot_exp;
}


/*Function is external*/
P= P + DALEC_LIKELIHOOD_GRACE_EWT(D);

/*EWT constraint*/
/*tot_exp=0;
double mewt=0, mewtm=0;
if (D.newt>0){*/
/*Relative volumetric constraint*/
/*Note: constraint imposed on t+1 of H2O pools*/
/*for (n=0;n<D.newt;n++){dn=D.ewtpts[n];mewtm+=(D.M_POOLS[D.nopools*dn+6] + D.M_POOLS[D.nopools*dn+7])/D.newt;mewt+=D.EWT[dn]/D.newt;}

for (n=0;n<D.newt;n++){dn=D.ewtpts[n];tot_exp+=pow((D.M_POOLS[D.nopools*dn+6]+ D.M_POOLS[D.nopools*dn+7]-D.EWT[dn]-mewtm+mewt)/D.ewt_obs_unc,2);}
P=P-0.5*tot_exp;}*/
/*
printf("mewt = %2.2f\n",mewt);
printf("mewtm = %2.2f\n",mewtm);
printf("totexp EWT = %2.2f\n",tot_exp);
*/


/*Note: only use with model ID = 806*/

if (D.ID==806){P = P + DALEC_806_MFCF(D);}
if (D.ID==807){P = P + DALEC_807_MFCF(D);}
if (D.ID==808){P = P + DALEC_807_MFCF(D);}







/*Constrain fire emissions here*/
if (D.otherpriors[2]>-9999){
/*Step 1. Sum fire emissions*/
for (n=0;n<D.nodays;n++){mfire+=D.M_FLUXES[n*D.nofluxes+16];}
/*Normalize fire constraint to daily mean flux*/
mfire=mfire/((double)D.nodays);

/*Step 2. Constrain against fire emissions*/
if (D.otherpriorunc[2]>0){P=P-0.5*pow(log(mfire/D.otherpriors[2])/log(D.otherpriorunc[2]),2);}
else {P=P-0.5*pow((mfire-D.otherpriors[2])/D.otherpriorunc[2],2);}
}



/*Constrain mean GPP*/
if (D.otherpriors[5]>-9999){mgpp=0;
/*Step 1. Sum fire emissions*/
for (n=0;n<D.nodays;n++){mgpp+=D.M_FLUXES[n*D.nofluxes+0];}
/*Normalize gpp constraint to daily mean flux*/
mgpp=mgpp/((double)D.nodays);
/*Step 2. Constrain against gpp*/
if (D.otherpriorunc[5]>0){P=P-0.5*pow(log(mgpp/D.otherpriors[5])/log(D.otherpriorunc[5]),2);}
else {P=P-0.5*pow((mgpp-D.otherpriors[5])/D.otherpriorunc[5],2);}
/*P=P-0.5*pow((mgpp-D.otherpriors[5])/D.otherpriorunc[5],2);*/

}


/*printf("MGPPprob = %2.2f\n",P);&*/



/*Constrain CMS disturbance fluxes*/
if (D.otherpriors[6]>-9999){mnpdf=0;
/*Step 1. Sum of biomass -> litter emissions (including fire)*/
for (n=0;n<D.nodays;n++){
for (f=0;f<7;f++){mnpdf+=D.M_FLUXES[n*D.nofluxes+npdfi[f]];}}
/*Normalize npdf constraint to daily mean flux*/
mnpdf=mnpdf/((double)D.nodays);
/*Step 2. Constrain against npdf*/
/*Only constrain if model flux < observed flux (otherwise mortality can explain internal cardamom fluxes)*/
if (mnpdf<D.otherpriors[6]) {
P=P-0.5*pow(log(mnpdf/D.otherpriors[6])/log(D.otherpriorunc[6]),2);
}}






if (isnan(P)){P=log(0);}
return P;}

