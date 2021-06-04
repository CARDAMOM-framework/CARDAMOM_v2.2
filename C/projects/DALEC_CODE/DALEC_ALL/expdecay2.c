#pragma once


/*Component of expdecay (next function)*/
double expdecay2(double const *POOLS, int pool, int nodays,int deltat,int nopools)
/*only accepting pool and number of days, the rest done here*/
{
/*using 365 day averaging window for each year!*/

/*explicitly calculating number of years*/
int noyears=floor(nodays*deltat/365);
int EDC=1,n,count,y;

/*yearly means and yearly means with offset Ds = 100 days*/
double P0=POOLS[pool];
/*Deriving exponential decay coefficients a,b and c in equation
 * Cexp = a + b*exp(c*t)*/
/*four mean pool values to be derived*/
/* MP0 = mean pool (year 1 to year end-2)
 * MP1 = mean pool (year 2 to year end-1)
 * MP0os = mean pool (year 1+os to year end-2+os)
 * MP1os = mean pool (year 1+os to year end-2+os)*/
double MP0=0,MP1=0,MP0os=0,MP1os=0;
/*averaging window*/
int aw=floor(365./(double)deltat);
/*offset in days*/
/*OFFSET = 1 is ideal to capture rapid exponential decay without compromising quality of fit*/
int os=1;
/*deriving coefficients to *
 * in Cexp = A + B exp(C*t);*/
double A,B,b,C;

/*mean pool within each averaging window*/
for (n=0;n<aw;n++){MP0=MP0+POOLS[n*nopools+pool];};MP0=MP0/(double)aw;
for (n=aw;n<aw*2;n++){MP1=MP1+POOLS[n*nopools+pool];};MP1=MP1/(double)aw;
for (n=os;n<aw+os;n++){MP0os=MP0os+POOLS[n*nopools+pool];};MP0os=MP0os/(double)aw;
for (n=aw+os;n<aw*2+os;n++){MP1os=MP1os+POOLS[n*nopools+pool];};MP1os=MP1os/(double)aw;

/*deriving mean gradient ratio dcdt1/dcdt0*/
/*dcdt1 is the numeric gradient between n+1 and n+365+1*/
/*dcdt0 is the numeric gradient between n and n+365*/
double dcdtr=0,dcdt1,dcdt0;
double dcdty1=(MP1os-MP0os);
double dcdty0=(MP1-MP0);
/*denominators*/

/*using multiyear mean to determine c*/
C=log(dcdty1/dcdty0)/((double)os*deltat);
/*deriving final exp. decay fit with startpoint = startpoint of pool*/

/*compared and validated against dalec_expdecay3.m*/
/*
printf("noyears = %d\n",noyears);
printf("MP0 = %f\n",MP0);
printf("MP1 = %f\n",MP1);
printf("MP0os = %f\n",MP0os);
printf("MP1os = %f\n",MP1os);
printf("A = %e\n",A);
printf("B = %e\n",B);
printf("C = %e\n",C);
printf("aw = %d\n",aw);
printf("dcdty1 = %f\n",dcdt1);
printf("dcdty0 = %f\n",dcdt0);
*/
/*half life must be more than noyears*/
/*if (fabs(-log(2)/C)<noyears*365.25 & C<0 & finite(C) & abs(B)>abs(A)*0.01){EDC=0;}*/
return C;}

