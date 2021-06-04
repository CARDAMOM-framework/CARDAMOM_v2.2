#include "CARDAMOM_DATA_STRUCTURE.c"
int ADD_PARAMETER_TO_STACK(DATA * DATA, char *parname,double parmin,double parmax){
/*TO DO: add "parnames" to DATA*/
/*Use realloc*/
/*Step 1. add more parameters*/
int N=DATA->nopars;

/*Step 2. re-allocate memory to accomodate more parameters*/

DATA->parname=realloc(DATA->parname,(N+1)*sizeof(char*));
DATA->parmin=realloc(DATA->parmin,(N+1)*sizeof(double));
DATA->parmax=realloc(DATA->parmax,(N+1)*sizeof(double));

/*Populate with new parameter entry*/
DATA->parmin[N]=parmin;
DATA->parmax[N]=parmax;
DATA->parname[N]=parname;
DATA->nopars=N+1;

return 0;}
