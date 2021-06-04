#pragma once
#include "../../../math_fun/std.c"
#include "../../../math_fun/mean.c"
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
for (n=0;n<50;n++){if (DATA.parpriors[n]>-9999 & DATA.parpriorunc[n]>-9999){p=p-0.5*pow(log(PARS[n]/DATA.parpriors[n])/log(DATA.parpriorunc[n]),2);}}

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

double mam=0;
double am=0;



if (D.ngpp>0){
tot_exp=0;
/*gppabs option = 1; directly assimilating GPP*/
/*gppabs option = 0; assimilating relative GPP constraint (e.g. SIF)*/
if (D.gppabs==1){
/*GPP*/
for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];tot_exp+=pow((D.M_GPP[dn]-D.GPP[dn])/2,2);}
P=P-0.5*tot_exp;}
else if (D.gppabs==0){
/*GPP relative*/
/*Step 1. Calculate mean at available datapoints*/
for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];mgpp+=D.M_GPP[dn];msif+=D.GPP[dn];
/*printf("CGPP(%i) = %4.4f\n", n+1, D.M_GPP[dn]); */
}
mgpp=mgpp/(double)D.ngpp;
msif=msif/(double)D.ngpp;
for (n=0;n<D.ngpp;n++){
dn=D.gpppts[n];tot_exp+=pow(log(D.M_GPP[dn]/mgpp/(D.GPP[dn]/msif))/log(2),2);


}
P=P-0.5*tot_exp;
}}

/*printf("GPPprob = %2.2f\n",P);*/


/*LAI likelyhood*/
tot_exp=0;
if (D.nlai>0 & D.otherpriors[4]<0){for (n=0;n<D.nlai;n++){dn=D.laipts[n];tot_exp+=pow(log(D.M_LAI[dn]/D.LAI[dn])/log(2),2);}
P=P-0.5*tot_exp;}
/*use timesteps for mean LAI calculation*/
else if (D.nlai>0 & D.otherpriors[4]>0){
for (n=0;n<D.nlai;n++){dn=D.laipts[n];mlai+=D.M_LAI[dn];}
P=P-0.5*pow((mlai/(double)D.nlai-D.otherpriors[4])/D.otherpriorunc[4],2);}
else if (D.otherpriors[4]>0){for (n=0;n<D.nodays;n++){mlai+=D.M_LAI[n];};mlai=mlai/(double)D.nodays;
/*P=P-0.5*pow((mlai-D.otherpriors[4])/D.otherpriorunc[4],2);}*/
P=P-0.5*pow(log(mlai/D.otherpriors[4])/log(D.otherpriorunc[4]),2);}
/*
printf("mlai = %2.2f\n",mlai);
printf("prob = %2.2f\n",-0.5*pow((mlai-D.otherpriors[4])/D.otherpriorunc[4],2));
*/

/*printf("LAIprob = %2.2f\n",P);*/


/*ET likelyhood*/
if (D.net>0){
tot_exp=0;
if (D.etiav<1){
/*Standard NEE likelihood approach*/
for (n=0;n<D.net;n++){dn=D.etpts[n];tot_exp+=pow(log(D.M_FLUXES[dn*D.nofluxes+28]/D.ET[dn])/log(2),2);}}
else{
/*Decoupling seasonal from interannual variations*/
/*Only use with monthly resolution fluxes, complete years & no missing data*/
/*Step 1. Mean model & data annual NBE*/
/*Step 2. Compare means*/
/*Step 3. Remove means from months for cost function*/
for (m=0;m<D.net/12;m++){
/*Calculate annual mean*/
mam=0;am=0;
for (n=0;n<12;n++){dn=D.etpts[n+m*12];mam=mam+D.M_FLUXES[dn*D.nofluxes+28];am=am+D.ET[dn];}
/*Calculate seasonal cost function*/
for (n=0;n<12;n++){dn=D.etpts[n+m*12];tot_exp+=pow(log(D.M_FLUXES[dn*D.nofluxes+28]/D.ET[dn]*am/mam)/log(2),2);}
/*Calculate annual cost function*/
tot_exp+=pow(log(am/mam)/log(1.2),2);
}}
P=P-0.5*tot_exp;}




/*printf("ETprob = %2.2f\n",P);*/




/*NEE likelyhood*/
if (D.nnee>0){
tot_exp=0;
if (D.neeiav<0){
/*Standard NEE likelihood approach*/
for (n=0;n<D.nnee;n++){dn=D.neepts[n];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn])/D.nee_obs_unc,2);}}
else{
/*Decoupling seasonal from interannual variations*/
/*Only use with monthly resolution fluxes, complete years & no missing data*/
/*Step 1. Mean model & data annual NBE*/
/*Step 2. Compare means*/
/*Step 3. Remove means from months for cost function*/
for (m=0;m<D.nnee/12;m++){
/*Calculate annual mean*/
mam=0;am=0;
for (n=0;n<12;n++){dn=D.neepts[n+m*12];mam=mam+D.M_NEE[dn];am=am+D.NEE[dn];}
/*normalize means*/
mam=mam/12;am=am/12;
/*Calculate seasonal cost function*/
for (n=0;n<12;n++){dn=D.neepts[n+m*12];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn]-mam+am)/D.nee_obs_unc,2);}
/*Calculate annual cost function*/
/*TEST: normalize model likelihood by normal distribution with mean zero and unc = x2 annual unc.*/
tot_exp+=pow((am-mam)/D.nee_annual_unc,2);


}}
P=P-0.5*tot_exp;}


/*printf("NEEprob = %2.2f\n",P);*/

/*Cwood likelyhood*/
/*continue from HERE!!!!*/
/*fraction of woody pool (with uncertainty of 2)*/
/*Change made on Jun 21 2016*/
/*for (n=1;n<D.nwoo;n++){dn=D.woopts[n];tot_exp+=pow(log((D.M_POOLS[D.nopools*dn+3]/D.M_POOLS[D.nopools*D.woopts[0]+3])/D.WOO[dn])/((D.WOO[dn]-D.WOO[D.woopts[0]])*log(2)),2);}*/
tot_exp=0;
if (D.nwoo>0){
for (n=0;n<D.nwoo;n++){dn=D.woopts[n];tot_exp+=pow(log(D.M_POOLS[D.nopools*dn+3]/D.WOO[dn])/log(D.otherpriorunc[1]),2);}
P=P-0.5*tot_exp;}


/*EWT constraint*/
tot_exp=0;
double mewt=0, mewtm=0;
if (D.newt>0){
/*Relative volumetric constraint*/
/*Note: constraint imposed on t+1 of H2O pools*/
for (n=0;n<D.newt;n++){dn=D.ewtpts[n];mewtm+=(D.M_POOLS[D.nopools*dn+6] + D.M_POOLS[D.nopools*dn+7])/D.newt;mewt+=D.EWT[dn]/D.newt;}

for (n=0;n<D.newt;n++){dn=D.ewtpts[n];tot_exp+=pow((D.M_POOLS[D.nopools*dn+6]+ D.M_POOLS[D.nopools*dn+7]-D.EWT[dn]-mewtm+mewt)/50,2);}
P=P-0.5*tot_exp;}
/*
printf("mewt = %2.2f\n",mewt);
printf("mewtm = %2.2f\n",mewtm);
printf("totexp EWT = %2.2f\n",tot_exp);
*/


/*Note: only use with model ID = 806*/

tot_exp=0;
double mband=0, mbandm=0,mwoo=0,mfol=0,mobs=0;
double modband=0, smod=0, sobs=0, amod=0, aobs=0;
double *bandmod;


if (D.ID==806){
if (D.nband1>0){
/*Constraint on relative magnitudes of foliar and wood pools*/
/*Step 1. Normalize time series*/
/*Step 2. weigh and add time-series (MOD)*/
/*Step 3. Compare normalized MOD anomalies vs normalized OBS anomalies*/
bandmod=calloc(D.nband1,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
mfol+=D.M_POOLS[D.nopools*dn + 1];
mwoo+=D.M_POOLS[D.nopools*dn + 3];
mobs+=D.BAND1[dn];}
mfol=mfol/(double)D.nband1;
mwoo=mwoo/(double)D.nband1;
mobs=mobs/(double)D.nband1;

/*weighed sum of timeseries*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
bandmod[n]=D.M_POOLS[D.nopools*dn + 1]*D.M_PARS[33]/mfol + D.M_POOLS[D.nopools*dn + 3]*(1-D.M_PARS[33])/mwoo;}

/*Derive MOD and OBS standard deviations*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
smod+=pow(bandmod[n]-1,2);
sobs+=pow(D.BAND1[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband1 -1));
sobs = sqrt(sobs/(double)(D.nband1 -1));

/*Derive normalized residuals*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
/*calculate standard anomalies*/
aobs=(D.BAND1[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}

if (D.nband2>0){
/*Constraint on relative magnitudes of foliar and wood pools*/
/*Step 1. Normalize time series*/
/*Step 2. weigh and add time-series (MOD)*/
/*Step 3. Compare normalized MOD anomalies vs normalized OBS anomalies*/
bandmod=calloc(D.nband2,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
mfol+=D.M_POOLS[D.nopools*dn + 1];
mwoo+=D.M_POOLS[D.nopools*dn + 3];
mobs+=D.BAND2[dn];}
mfol=mfol/(double)D.nband2;
mwoo=mwoo/(double)D.nband2;
mobs=mobs/(double)D.nband2;

/*weighed sum of timeseries*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
bandmod[n]=D.M_POOLS[D.nopools*dn + 1]*D.M_PARS[34]/mfol + D.M_POOLS[D.nopools*dn + 3]*(1-D.M_PARS[34])/mwoo;}

/*Derive MOD and OBS standard deviations*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
smod+=pow(bandmod[n]-1,2);
sobs+=pow(D.BAND2[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband2 -1));
sobs = sqrt(sobs/(double)(D.nband2 -1));

/*Derive normalized residuals*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
/*calculate standard anomalies*/
aobs=(D.BAND2[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}
}






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
P=P-0.5*pow(log(mgpp/D.otherpriors[5])/log(D.otherpriorunc[5]),2);
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

