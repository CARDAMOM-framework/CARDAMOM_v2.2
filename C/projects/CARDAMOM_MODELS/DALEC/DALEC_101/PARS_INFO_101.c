/*PARAMETER_INFO (typedef struct) must have at least 3 fields
 *  * npars,
 *   * parmax
 *    * parmin*/
/*where is it defined?*/
/*For DALEC_FIRES: as all GPP allocation fractions are inter-dependent (sum = 1)*/
/*MCMC sampling of GPP allocation priors approximated as 0.01-0.5 NPP for*/
/*photosynthetic pools and 0.01-1 of remaining NPP for root and wood pool*/

int PARS_INFO_101(DATA *CARDADATA)
{

CARDADATA->nopars=14;

/*Total ecosystem RT*/
CARDADATA->parmin[0]=0.000001;
CARDADATA->parmax[0]=0.1;

/*Fraction of GPP respired*/
CARDADATA->parmin[1]=0.2;
CARDADATA->parmax[1]=0.8;

/*Leaf Lifespan*/
CARDADATA->parmin[2]=1.001;
CARDADATA->parmax[2]=8;

/*Temp factor* = Q10 = 1.2(0.018)-1.6(0.08)*/
CARDADATA->parmin[3]=0.018;
CARDADATA->parmax[3]=0.08;

/*Canopy Efficiency*/
CARDADATA->parmin[4]=5;
CARDADATA->parmax[4]=50;

/*Bday*/
CARDADATA->parmin[5]=365.25;
CARDADATA->parmax[5]=365.25*4;

/*Fraction to Clab*/
CARDADATA->parmin[6]=0.01;
CARDADATA->parmax[6]=0.5;

/*Clab Release period*/
CARDADATA->parmin[7]=10;
CARDADATA->parmax[7]=100;

/*Fday*/
CARDADATA->parmin[8]=365.25;
CARDADATA->parmax[8]=365.25*4;

/*Leaf fall period*/
CARDADATA->parmin[9]=20;
CARDADATA->parmax[9]=150;

/*LMCA*/
/*Kattge et al. 2011*/
/*Kattge et al., provide a range of 10 400 g m-2, i.e. 5 200 gC m-2*/
CARDADATA->parmin[10]=5;
CARDADATA->parmax[10]=200;

/*INITIAL VALUES DECLARED HERE*/

/*C labile*/
CARDADATA->parmin[11]=1.0;
CARDADATA->parmax[11]=2000.0;

/*C foliar*/
CARDADATA->parmin[12]=1.0;
CARDADATA->parmax[12]=2000.0;

/*C other*/
CARDADATA->parmin[13]=1.0;
CARDADATA->parmax[13]=200000.0;

return 0;

}


