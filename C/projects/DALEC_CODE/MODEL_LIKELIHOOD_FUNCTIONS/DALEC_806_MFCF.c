/*Multi-frequency cost function*/
double DALEC_806_MFCF(DATA D){

double mband=0, mbandm=0,mwoo=0,mfol=0,mobs=0;
double modband=0, smod=0, sobs=0, amod=0, aobs=0;
double *bandmod, P=1, tot_exp;
int n, dn;

if (D.nband1>0){
tot_exp=0;
/*Constraint on relative magnitudes of foliar and wood pools*/
/*Step 1. Normalize time series*/
/*Step 2. weigh and add time-series (MOD)*/
/*Step 3. Compare normalized MOD anomalies vs normalized OBS anomalies*/
bandmod=calloc(D.nband1,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
/*M_POOLS includes initial conditions, offsetting time index by +1*/
mfol+=D.M_POOLS[D.nopools*(dn+1) + 1];
mwoo+=D.M_POOLS[D.nopools*(dn+1) + 3];
mobs+=D.BAND1[dn];}
mfol=mfol/(double)D.nband1;
mwoo=mwoo/(double)D.nband1;
mobs=mobs/(double)D.nband1;

/*weighed sum of timeseries*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
bandmod[n]=D.M_POOLS[D.nopools*(dn+1) + 1]*D.M_PARS[33]/mfol + D.M_POOLS[D.nopools*(dn+1) + 3]*(1-D.M_PARS[33])/mwoo;}

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
tot_exp=0;
/*Constraint on relative magnitudes of foliar and wood pools*/
/*Step 1. Normalize time series*/
/*Step 2. weigh and add time-series (MOD)*/
/*Step 3. Compare normalized MOD anomalies vs normalized OBS anomalies*/
bandmod=calloc(D.nband2,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
mfol+=D.M_POOLS[D.nopools*(dn+1) + 1];
mwoo+=D.M_POOLS[D.nopools*(dn+1) + 3];
mobs+=D.BAND2[dn];}
mfol=mfol/(double)D.nband2;
mwoo=mwoo/(double)D.nband2;
mobs=mobs/(double)D.nband2;

/*weighed sum of timeseries*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
bandmod[n]=D.M_POOLS[D.nopools*(dn+1) + 1]*D.M_PARS[34]/mfol + D.M_POOLS[D.nopools*(dn+1) + 3]*(1-D.M_PARS[34])/mwoo;}

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


return P;}





