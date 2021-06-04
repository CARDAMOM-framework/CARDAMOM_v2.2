double DALEC_LIKELIHOOD_NEE(DATA D){



int dn, n,m;
double am,mam,P=0;
/*NEE likelyhood*/
if (D.nnee>0){
double tot_exp=0;
if (D.nee_annual_unc<0){

if (D.nneeunc==0){
/*Standard NEE likelihood approach*/
for (n=0;n<D.nnee;n++){dn=D.neepts[n];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn])/D.nee_obs_unc,2);}}


if (D.nneeunc>0){
/*Standard NEE likelihood approach*/
for (n=0;n<D.nnee;n++){dn=D.neepts[n];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn])/D.NEEunc[dn],2);}}

}
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



return P;
}
