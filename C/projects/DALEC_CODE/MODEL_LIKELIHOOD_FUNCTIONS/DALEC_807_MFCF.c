/*Multi-frequency cost function*/
double DALEC_807_MFCF(DATA D){

/*7 parameters*/
/*PARS[33 - 39]*/
/*Pband (band 1)*/
double Pband_h2o_rz=D.M_PARS[33];
double Pband_h2o_ns=(1 - Pband_h2o_rz)*D.M_PARS[34];
double Pband_woo=1-Pband_h2o_rz-Pband_h2o_ns;
/*Lband (band 2)*/
double Lband_h2o_rz=D.M_PARS[35];
double Lband_h2o_ns=(1 - Lband_h2o_rz)*D.M_PARS[36];
double Lband_woo=(1-Lband_h2o_rz-Lband_h2o_ns)*D.M_PARS[37];
double Lband_fol=1- Lband_h2o_rz - Lband_h2o_ns - Lband_woo;
/*Cband (band 3)*/
double Cband_h2o_ns=D.M_PARS[38];
double Cband_woo=(1-Cband_h2o_ns)*D.M_PARS[39];
double Cband_fol=1-Cband_woo - Cband_h2o_ns;
/*Kuband (band 4)*/
double Kuband_fol=D.M_PARS[40];
double Kuband_prec=1-Kuband_fol;



double mband, mbandm,mwoo,mfol,mobs,mprec,mh2o_rz,mh2o_ns;
double modband, smod, sobs, amod, aobs;
double *bandmod, P=0, tot_exp;
int n, dn;


/*number of MET drivers*/
int nomet=((DALEC *)D.MODEL)->nomet;

/*number of DALEC pools*/
int nopools=((DALEC *)D.MODEL)->nopools;



/*Constraint on relative magnitudes of foliar and wood pools*/
/*Step 1. Normalize time series*/
/*Step 2. weigh and add time-series (MOD)*/
/*Step 3. Compare normalized MOD anomalies vs normalized OBS anomalies*/


/*Pband*/
if (D.nband1>0){
mh2o_ns=0;mh2o_rz=0;mwoo=0;mobs=0;mfol=0;mprec=0;smod=0;sobs=0;tot_exp=0;
bandmod=calloc(D.nband1,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
/*M_POOLS includes initial conditions, offsetting time index by +1*/
mh2o_ns+=D.MET[nomet*dn + 8];
mh2o_rz+=D.M_POOLS[nopools*(dn+1) + 6];
mwoo+=D.M_POOLS[nopools*(dn+1) + 3] + D.M_POOLS[nopools*(dn+1) + 0];
mobs+=D.BAND1[dn];}
mh2o_ns=mh2o_ns/(double)D.nband1;
mh2o_rz=mh2o_rz/(double)D.nband1;
mwoo=mwoo/(double)D.nband1;
mobs=mobs/(double)D.nband1;
/*weighed sum of timeseries*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
bandmod[n]+=D.MET[nomet*dn + 8]*Pband_h2o_ns/mh2o_ns;
bandmod[n]+=D.M_POOLS[D.nopools*(dn+1) + 6]*Pband_h2o_rz/mh2o_rz;
bandmod[n]+=(D.M_POOLS[D.nopools*(dn+1) + 3] + D.M_POOLS[D.nopools*(dn+1) + 0])*Pband_woo/mwoo;}
/*Derive MOD and OBS standard deviations (same for each band)*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];smod+=pow(bandmod[n]-1,2);sobs+=pow(D.BAND1[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband1 -1));
sobs = sqrt(sobs/(double)(D.nband1 -1));

/*Derive normalized residuals (and calculate standard anomalies)*/
for (n=0;n<D.nband1;n++){dn=D.band1pts[n];
aobs=(D.BAND1[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}


/*Lband*/
if (D.nband2>0){
mh2o_ns=0;mh2o_rz=0;mwoo=0;mobs=0;mfol=0;mprec=0;smod=0;sobs=0;tot_exp=0;
bandmod=calloc(D.nband2,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
/*M_POOLS includes initial conditions, offsetting time index by +1*/
mh2o_ns+=D.MET[nomet*dn + 8];
mh2o_rz+=D.M_POOLS[nopools*(dn+1) + 6];
mwoo+=D.M_POOLS[nopools*(dn+1) + 3] + D.M_POOLS[nopools*(dn+1) + 0];
mfol+=D.M_POOLS[nopools*(dn+1) + 1];
mobs+=D.BAND2[dn];}
mh2o_ns=mh2o_ns/(double)D.nband2;
mh2o_rz=mh2o_rz/(double)D.nband2;
mwoo=mwoo/(double)D.nband2;
mfol=mfol/(double)D.nband2;
mobs=mobs/(double)D.nband2;
/*weighed sum of timeseries*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
bandmod[n]+=D.MET[nomet*dn + 8]*Lband_h2o_ns/mh2o_ns;
bandmod[n]+=D.M_POOLS[D.nopools*(dn+1) + 6]*Lband_h2o_rz/mh2o_rz;
bandmod[n]+=(D.M_POOLS[D.nopools*(dn+1) + 3] + D.M_POOLS[D.nopools*(dn+1) + 0])*Lband_woo/mwoo;
bandmod[n]+=D.M_POOLS[D.nopools*(dn+1) + 1]*Lband_fol/mfol;}
/*Derive MOD and OBS standard deviations (same for each band)*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];smod+=pow(bandmod[n]-1,2);sobs+=pow(D.BAND2[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband2 -1));
sobs = sqrt(sobs/(double)(D.nband2 -1));

/*Derive normalized residuals (and calculate standard anomalies)*/
for (n=0;n<D.nband2;n++){dn=D.band2pts[n];
aobs=(D.BAND2[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}


/*Cband*/
if (D.nband3>0){
mh2o_ns=0;mh2o_rz=0;mwoo=0;mobs=0;mfol=0;mprec=0;smod=0;sobs=0;tot_exp=0;
bandmod=calloc(D.nband3,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband3;n++){dn=D.band3pts[n];
/*M_POOLS includes initial conditions, offsetting time index by +1*/
mh2o_ns+=D.MET[nomet*dn + 8];
mwoo+=D.M_POOLS[nopools*(dn+1) + 3] + D.M_POOLS[nopools*(dn+1) + 0];
mfol+=D.M_POOLS[nopools*(dn+1) + 1];
mobs+=D.BAND3[dn];}
mh2o_ns=mh2o_ns/(double)D.nband3;
mwoo=mwoo/(double)D.nband3;
mfol=mfol/(double)D.nband3;
mobs=mobs/(double)D.nband3;
/*weighed sum of timeseries*/
for (n=0;n<D.nband3;n++){dn=D.band3pts[n];
bandmod[n]+=D.MET[nomet*dn + 8]*Cband_h2o_ns/mh2o_ns;
bandmod[n]+=(D.M_POOLS[D.nopools*(dn+1) + 3] + D.M_POOLS[D.nopools*(dn+1) + 0])*Cband_woo/mwoo;
bandmod[n]+=D.M_POOLS[D.nopools*(dn+1) + 1]*Cband_fol/mfol;}
/*Derive MOD and OBS standard deviations (same for each band)*/
for (n=0;n<D.nband3;n++){dn=D.band3pts[n];smod+=pow(bandmod[n]-1,2);sobs+=pow(D.BAND3[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband3 -1));
sobs = sqrt(sobs/(double)(D.nband3 -1));
/*Derive normalized residuals (and calculate standard anomalies)*/
for (n=0;n<D.nband3;n++){dn=D.band3pts[n];
aobs=(D.BAND3[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}



/*Kuband*/
if (D.nband4>0){
mh2o_ns=0;mh2o_rz=0;mwoo=0;mobs=0;mfol=0;mprec=0;smod=0;sobs=0;tot_exp=0;
bandmod=calloc(D.nband4,sizeof(double));
/*Derive means*/
for (n=0;n<D.nband4;n++){dn=D.band4pts[n];
/*M_POOLS includes initial conditions, offsetting time index by +1*/
mprec+=D.MET[nomet*dn + 8];
mfol+=D.M_POOLS[nopools*(dn+1) + 1];
mobs+=D.BAND4[dn];}
mprec=mprec/(double)D.nband4;
mfol=mfol/(double)D.nband4;
mobs=mobs/(double)D.nband4;
/*weighed sum of timeseries*/
for (n=0;n<D.nband4;n++){dn=D.band4pts[n];
bandmod[n]+=D.MET[nomet*dn + 8]*Kuband_prec/mprec;
bandmod[n]+=D.M_POOLS[D.nopools*(dn+1) + 1]*Kuband_fol/mfol;}
/*Derive MOD and OBS standard deviations (same for each band)*/
for (n=0;n<D.nband4;n++){dn=D.band4pts[n];smod+=pow(bandmod[n]-1,2);sobs+=pow(D.BAND4[dn]-mobs,2);}
smod = sqrt(smod/(double)(D.nband4 -1));
sobs = sqrt(sobs/(double)(D.nband4 -1));
/*Derive normalized residuals (and calculate standard anomalies)*/
for (n=0;n<D.nband4;n++){dn=D.band4pts[n];
aobs=(D.BAND4[dn] - mobs)/sobs;
amod=(bandmod[n] - 1)/smod;
tot_exp+=pow((aobs - amod)/1,2);}
P=P-0.5*tot_exp;
free(bandmod);}










return P;}





