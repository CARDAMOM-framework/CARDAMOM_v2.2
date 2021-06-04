#pragma once
int EDC1_831(double const *pars, DATA DATA, struct EDCDIAGNOSTIC *EDCD)
{

/*EDCD=EDCD2;*/
/*List of EDCs in this file*/
/*1.  TOR_LIT faster than TOR_SOM
2. Litter2SOM greater than SOM to Atm. rate
3. TOR_FOL faster than TOR_WOOD
4. Root turnover greater than SOM turnover at meantemp
5. Allocation to canopy is comparable to allocation to fine roots
6-14. Taken up by EDC2
15. Foliage CF> wood CF & Foliage CF > Soil CF 
*/



double meantemp=DATA.meantemp;
double meanrad=DATA.meanrad;

/*This function was created on 7 Jan 2014*/
/*Rules here are specified as in Bloom et al., 2014 paper*/
/*EDC1 contains all checks that do not require a full DALEC_CDEA (DALEC2) model run*/


/*ALL EDCs set as 1 initially*/
EDCD->nedc=100;
int n; for (n=0;n<EDCD->nedc;n++){EDCD->PASSFAIL[n]=1;}
/*declaring variables and constants for subsequent EDCs*/
int EDC=1;
int DIAG=EDCD->DIAG;
/*obsolete edcc constant was kept here*/
/*deriving true allocation fractions*/
double const fauto=pars[1];
double const ffol=(1-fauto)*pars[2];
double const flab=(1-fauto-ffol)*pars[12];
double const froot=(1-fauto-ffol-flab)*pars[3];
double const fwood=1-fauto-ffol-flab-froot;
/*fraction of GPP som under equilibrium conditions*/
double const fsom=fwood+(froot+flab+ffol)*pars[0]/(pars[0]+pars[7]);

/*yearly leaf loss fraction*/
double torfol=1/(pars[4]*365.25);


/*EDC CHECK NO 1*/
/*TOR_LIT faster than TOR_SOM*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[1-1]==1)) & (pars[8]>pars[7])){EDC=0;EDCD->PASSFAIL[1-1]=0;}

/*EDC CHECK NO 2*/
/*Litter2SOM greater than SOM to Atm. rate*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[2-1]==1)) & (pars[0]<pars[8])){EDC=0;EDCD->PASSFAIL[2-1]=0;}

/*EDC CHECK NO 3*/
/*TOR_FOL faster than TOR_WOOD */
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[3-1]==1)) & (pars[5]>torfol)){EDC=0;EDCD->PASSFAIL[3-1]=0;}

/*EDC CHECK NO 4*/
/*Root turnover greater than SOM turnover at meantemp*/
/*same as this*/
/*\text{EDC 4: }(1-\pavii)^{365} > \Pi_{i=1}^{365} (1-\paix \tratei)*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[4-1]==1)) & (pars[6]<pars[8]*exp(pars[9]*meantemp))){EDC=0;EDCD->PASSFAIL[4-1]=0;}

/*EDC no 5 is addressed in EDC2_FIREBUCKET.c*/

/*EDC CHECK NO 6*/
/*Bday Fday difference>45 days */
/*if (((EDC==1 & DIAG==0) || DIAG==1) & (fabs(pars[14] - pars[11])<45 | fabs(pars[14]-pars[11])>320.25)){EDC=0;EDCD->PASSFAIL[6-1]=0;}*/
/*now obsolete!*/

/*EDC CHECK NO 5*/
/*Allocation to canopy is comparable to allocation to fine roots*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[5-1]==1)) & ((ffol+flab)>5*froot | (ffol+flab)*5<froot)){EDC=0;EDCD->PASSFAIL[5-1]=0;}

/*EDC No 15* Foliage CF> wood CF & SOM CF*/
if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[15-1]==1)) & (pars[27]<pars[28] | pars[27]<pars[29])){EDC=0;EDCD->PASSFAIL[15-1]=0;}

/*EDC No 16: Field capacity > wilting point*/
/*if (((EDC==1 & DIAG==0) || DIAG==1 || (EDC==1 & DIAG==2 & EDCD->SWITCH[16-1]==1)) & (pars[24]<pars[25])){EDC=0;EDCD->PASSFAIL[16-1]=0;}
*/


/*Add any generalisations derivable from EDC2 (post-run checks) here*/
/*Note: these must be tested to ensure that DALEC2 run is NOT needed */






return EDC;


}






