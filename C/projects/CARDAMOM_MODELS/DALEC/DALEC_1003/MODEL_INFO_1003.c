#pragma once
#include "PARS_INFO_1003.c"
#include "DALEC_1003.c"
#include "EDC1_1003.c"
#include "EDC2_1003.c"
#include "../../../DALEC_CODE/MODEL_LIKELIHOOD_FUNCTIONS/DALEC_MLF.c"


int MODEL_INFO_1003(DATA * DATA){

/*Step 1. Declare structure*/
/*"static" ensures that the memory is declared in one instance and visible to all functions (I think)*/
static DALEC DALECmodel;

/*Step 2: Fill structure with model-specific info*/
DALECmodel.nopools=8;
DALECmodel.nomet=9;/*This should be compatible with CBF file, if not then disp error*/
DALECmodel.nopars=37;
DALECmodel.nofluxes=33;

/*Short-term: copy quantities into DATA structure to reduce dependencies in CARDAMOM_MODEL_LIBRARY.c*/
/*Long-term: remove dependencies on DATA.nofluxes... etc. in CARDAMOM_READ_BINARY_DATA and DALEC_ALL_LIKELIHOOD.c*/
DATA->nopools=DALECmodel.nopools;
DATA->nopars=DALECmodel.nopars;
DATA->nofluxes=DALECmodel.nofluxes;

/*All model functions*/
/*User is able to add further functions as deemed necessary*/
/*Function names are declared in ../DALEC_ALL/DALEC_MODULE.c*/
/*Consider starting new module for radically different model structures*/
DALECmodel.dalec=DALEC_1003;
DALECmodel.edc1=EDC1_1003;
DALECmodel.edc2=EDC2_1003;

/*Initialize parameter fields*/
/*initializing parmin and parmax fields*/
/*Currently assigned to "DATA", since MCMC needs info separately*/
INITIALIZE_PARAMETER_FIELDS(DATA);
PARS_INFO_1003(DATA);

oksofar("about to declare EDCD");
printf("DALECmodel.EDCD = %p\n",DALECmodel.EDCD);
/*Initialize the EDCD structure*/
EDCSETUP(*DATA,&DALECmodel.EDCD);
oksofar("done with declaration");
printf("DALECmodel.EDCD->EQF = %2.2f\n",DALECmodel.EDCD->EQF);



/*initializing model*/
DATA->MODEL=&DALECmodel;
DATA->MLF=DALEC_MLF;




return 0;}
