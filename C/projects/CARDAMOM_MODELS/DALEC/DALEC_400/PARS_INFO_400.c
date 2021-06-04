/*PARAMETER_INFO (typedef struct) must have at least 3 fields
 *  * npars,
 *   * parmax
 *    * parmin*/
/*where is it defined?*/
/*For DALEC_FIRES: as all GPP allocation fractions are inter-dependent (sum = 1)*/
/*MCMC sampling of GPP allocation priors approximated as 0.01-0.5 NPP for*/
/*photosynthetic pools and 0.01-1 of remaining NPP for root and wood pool*/

int PARS_INFO_400(DATA *CARDADATA)
{

CARDADATA->nopars=23;

/*Decomposition rate*/
CARDADATA->parmin[0]=0.00001;
CARDADATA->parmax[0]=0.01;

/*Fraction of GPP respired*/
CARDADATA->parmin[1]=0.2;
CARDADATA->parmax[1]=0.8;

/*Fraction to foliage*/
CARDADATA->parmin[2]=0.01;
CARDADATA->parmax[2]=0.5;

/*Fraction to roots*/
CARDADATA->parmin[3]=0.01;
CARDADATA->parmax[3]=1;

/*Leaf Lifespan*/
CARDADATA->parmin[4]=1.001;
CARDADATA->parmax[4]=8;


/*TOR wood*/
CARDADATA->parmin[5]=0.000025;
CARDADATA->parmax[5]=0.001;


/*TOR roots*/
CARDADATA->parmin[6]=0.0001;
CARDADATA->parmax[6]=0.01;

/*TOR litter*/
CARDADATA->parmin[7]=0.0001;
CARDADATA->parmax[7]=0.01;

/*TOR SOM*/
CARDADATA->parmin[8]=0.0000001;
CARDADATA->parmax[8]=0.001;


/*Temp factor*/
CARDADATA->parmin[9]=0.018;
CARDADATA->parmax[9]=0.08;


/*Canopy Efficiency*/
CARDADATA->parmin[10]=5;
CARDADATA->parmax[10]=50;

/*Bday*/
CARDADATA->parmin[11]=365.25;
CARDADATA->parmax[11]=365.25*4;

/*Fraction to Clab*/
CARDADATA->parmin[12]=0.01;
CARDADATA->parmax[12]=0.5;

/*Clab Release period*/
CARDADATA->parmin[13]=10;
CARDADATA->parmax[13]=100;

/*Fday*/
CARDADATA->parmin[14]=365.25;
CARDADATA->parmax[14]=365.25*4;

/*Leaf fall period*/
CARDADATA->parmin[15]=20;
CARDADATA->parmax[15]=150;

/*LMCA*/
CARDADATA->parmin[16]=5;
CARDADATA->parmax[16]=200;
/*Kattge et al. 2011*/
/*Kattge et al., provide a range of 10 400 g m-2, i.e. 5 200 gC m-2*/

/*INITIAL VALUES DECLARED HERE*/

/*C labile @ t=0*/
CARDADATA->parmin[17]=1.0;
CARDADATA->parmax[17]=2000.0;

/*C foliar @ t=0*/
CARDADATA->parmin[18]=1.0;
CARDADATA->parmax[18]=2000.0;

/*C roots @ t=0*/
CARDADATA->parmin[19]=1.0;
CARDADATA->parmax[19]=2000.0;

/*C wood @ t=0*/
CARDADATA->parmin[20]=1.0;
CARDADATA->parmax[20]=100000.0;

/*C litter @ t=0*/
CARDADATA->parmin[21]=1.0;
CARDADATA->parmax[21]=2000.0;

/*C som @ t=0*/
CARDADATA->parmin[22]=1.0;
CARDADATA->parmax[22]=200000.0;

return 0;

}


