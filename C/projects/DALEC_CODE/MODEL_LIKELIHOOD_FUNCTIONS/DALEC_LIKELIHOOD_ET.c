double DALEC_LIKELIHOOD_ET(DATA D){
/*Data structure, includes model and data*/

/*EWT constraint*/
double tot_exp=0;
double am=0, mam=0;
int n,m,dn;
double P=0;
/*General notes*/
/* If D.et_annual_unc<1, then ET constraint is occurring on monthly basis*/
/* For log_et_obs*/
double th=D.et_obs_threshold;

/*ET likelyhood*/
if (D.net>0){
tot_exp=0;
if (D.et_annual_unc<1){
/*Standard NEE likelihood approach*/
for (n=0;n<D.net;n++){dn=D.etpts[n];tot_exp+=pow(log(max(D.M_FLUXES[dn*D.nofluxes+28],th)/max(D.ET[dn],th))/log(D.et_obs_unc),2);}}
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
for (n=0;n<12;n++){dn=D.etpts[n+m*12];tot_exp+=pow(log(max(D.M_FLUXES[dn*D.nofluxes+28],th)/max(D.ET[dn],th)*am/mam)/log(D.et_obs_unc),2);}
/*Calculate annual cost function*/
tot_exp+=pow(log(am/mam)/log(D.et_annual_unc),2);
}}
P=P-0.5*tot_exp;}





return P;
}
