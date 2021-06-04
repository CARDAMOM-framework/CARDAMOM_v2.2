double DALEC_LIKELIHOOD_CH4(DATA D){
/*Data structure, includes model and data*/

/*EWT constraint*/
double tot_exp=0;
double am=0, mam=0;
int n,m,dn;
double P=0;
/*General notes*/
/* If D.et_annual_unc<1, then ET constraint is occurring on monthly basis*/
/* For log_et_obs*/
double th=D.ch4_obs_threshold;

/*CH4 likelyhood*/
if (D.nch4>0){
tot_exp=0;

if (D.ch4_annual_unc<1){
/*Standard NEE likelihood approach*/
for (n=0;n<D.nch4;n++){dn=D.ch4pts[n];tot_exp+=pow(log(max(D.M_FLUXES[dn*D.nofluxes+37],th)/max(D.CH4[dn],th))/log(D.ch4_obs_unc),2);}}
else{
/*Decoupling seasonal from interannual variations*/
/*Only use with monthly resolution fluxes, complete years & no missing data*/
/*Step 1. Mean model & data annual NBE*/
/*Step 2. Compare means*/
/*Step 3. Remove means from months for cost function*/
for (m=0;m<D.nch4/12;m++){
/*Calculate annual mean*/
mam=0;am=0;
for (n=0;n<12;n++){dn=D.ch4pts[n+m*12];mam=mam+D.M_FLUXES[dn*D.nofluxes+37];am=am+D.CH4[dn];}
/*Calculate seasonal cost function*/
for (n=0;n<12;n++){dn=D.ch4pts[n+m*12];tot_exp+=pow(log(max(D.M_FLUXES[dn*D.nofluxes+37],th)/max(D.CH4[dn],th)*am/mam)/log(D.ch4_obs_unc),2);}
/*Calculate annual cost function*/
tot_exp+=pow(log(am/mam)/log(D.ch4_annual_unc),2);
}}
P=P-0.5*tot_exp;}





return P;
}
