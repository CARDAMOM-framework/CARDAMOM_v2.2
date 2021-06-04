#pragma once
/*mean matrix from double pointer routine*/
double mean_annual_pool(double *POOLS, int year, int pool, int nopools,int deltat){
/*inputs
 * POOLS: Pools double pointer, as output from DALEC
 * year: year for which to average (first year = 0)
 * pool: the specific pool 
 * nc
declarations*/
int r,c;
double meanpool=0;
/*deriving mean of pool p*/
int stday=floor(365.25*year/deltat);
int enday=floor(365.25*(year+1)/deltat);
for (c=stday;c<enday;c++){
meanpool=meanpool+POOLS[c*nopools+pool]/(enday-stday);}
/*returing meanpool value*/
return meanpool;}
