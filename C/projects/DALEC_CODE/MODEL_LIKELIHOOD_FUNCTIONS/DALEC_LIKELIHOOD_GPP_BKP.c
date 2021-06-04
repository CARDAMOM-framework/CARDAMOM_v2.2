double DALEC_LIKELIHOOD_GPP(DATA D){
/*Data structure, includes model and data*/

/*EWT constraint*/
double tot_exp=0;
double am=0, mam=0;
int n,m,dn;
double P=0;
double mgpp=1,msif=1;
/*General notes*/
/* If D.et_annual_unc<1, then ET constraint is occurring on monthly basis*/
/* For log_et_obs*/
double thresh=D.gpp_obs_threshold;

/*GPP likelyhood*/
if (D.ngpp>0){
tot_exp=0;

/*Normalization option for SIF*/
if (D.gppabs==0){
mgpp=0,msif=0;
for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];mgpp+=D.M_GPP[dn];msif+=D.GPP[dn];}
mgpp=mgpp/(double)D.ngpp;
msif=msif/(double)D.ngpp;}


if (D.gpp_annual_unc<1){
/*Standard NEE likelihood approach*/
for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];tot_exp+=pow(log(D.M_GPP[dn]/mgpp/(D.GPP[dn]/msif))/log(max(1,thresh/(D.GPP[dn]/msif))*D.gpp_obs_unc),2);}}
else{
/*Decoupling seasonal from interannual variations*/
/*Only use with monthly resolution fluxes, complete years & no missing data*/
/*Step 1. Mean model & data annual NBE*/
/*Step 2. Compare means*/
/*Step 3. Remove means from months for cost function*/
for (m=0;m<D.ngpp/12;m++){
/*Calculate annual mean*/
mam=0;am=0;
for (n=0;n<12;n++){dn=D.gpppts[n+m*12];mam=mam+D.M_GPP[dn]/mgpp;am=am+D.GPP[dn]/msif;}
/*Calculate seasonal cost function*/
for (n=0;n<12;n++){dn=D.gpppts[n+m*12];tot_exp+=pow(log(D.M_GPP[dn]/mgpp/mam/(D.GPP[dn]/msif/am))/log(max(1,thresh/(D.GPP[dn]/msif))*D.gpp_obs_unc),2);}
/*Calculate annual cost function*/
tot_exp+=pow(log(am/mam)/log(D.gpp_annual_unc),2);
}}
P=P-0.5*tot_exp;}


printf("am = %2.2f\n",am);
printf("mam = %2.2f\n",mam);


return P;
}
