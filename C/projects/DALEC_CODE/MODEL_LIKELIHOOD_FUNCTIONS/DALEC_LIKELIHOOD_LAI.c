double DALEC_LIKELIHOOD_LAI(DATA D){


double tot_exp=0, P=0, mlai=0;
int n,dn;
/*LAI likelihood*/
/*For future reference: average of model LAI at t and t+1 likely more appropriate*/
/*
double *LAIA=calloc(D.nodays, sizeof(double));
*/
/*Step 2. calculate t<->t+1 LAI*/
/*for (n=0;n<D.nodays;n++){LAIA[n]=}
*/


/*D.LAI[dn] = observed LAI at index dn*/
/*D.M_LAI[dn] = modelled LAI at index dn*/

/*Standard timestep model-data comparison*/

//Notes: If no "MLAI" (mean LAI) is provided, then cost function conducts a month-to-month comparison
if (D.otherpriors[4]<0 & D.nlai>0){for (n=0;n<D.nlai;n++){dn=D.laipts[n];tot_exp+=pow(log(D.M_LAI[dn]/D.LAI[dn])/log(2),2);}
P=P-0.5*tot_exp;}

/*use timesteps for mean LAI calculation*/
/*Here only using mean LAI*/
//NOTES: if "MLAI" (mean LAI) is provided, then only mean LAI is considered (either throughout model run or throughout overlapping timeseries, depending on whether LAI timeseries data is provided.
else if (D.otherpriors[4]>0){

//If timeseries data is provided, calculate model MLAI constraint during those points
if (D.nlai>0){
for (n=0;n<D.nlai;n++){dn=D.laipts[n];mlai+=D.M_LAI[dn];}
mlai=mlai/(double)D.nlai;}
//If no timeseries data is provided, calculate model MLAI for whole run
else {for (n=0;n<D.nodays;n++){mlai+=D.M_LAI[n];};mlai=mlai/(double)D.nodays;}

//Apply MLAI as normal or log-normal distrubution (depending if MLAI unc is positive or negative)
if (D.otherpriorunc[4]<0){
P=P-0.5*pow((mlai-D.otherpriors[4])/D.otherpriorunc[4],2);}
else {
P=P-0.5*pow(log(mlai/D.otherpriors[4])/log(D.otherpriorunc[4]),2);}

}


return P;
}
