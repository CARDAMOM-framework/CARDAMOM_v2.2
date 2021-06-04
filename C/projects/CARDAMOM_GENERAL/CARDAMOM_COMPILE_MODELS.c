
#include "../../auxi_fun/okcheck.c"
#include "math.h"
#include <stdio.h>
#include "CARDAMOM_READ_BINARY_DATA.c"
/*syntax CARDAMOM_READ_BINARY_CARDADATA(char *filename,CARDADATA *CARDADATA)*/




/* the gateway function */
/*NOTE: this function can be used with any valid model ID (stored in CARDADATA.ID)*/
int main()
{
 
/*declaring loop variable n*/ 

/*storing command line inputs as 2 files*/



/*declaring data structure*/
DATA DATA;

  
CARDAMOM_MODEL_INFO(&DATA,1);

/*
CARDAMOM_READ_BINARY_DATA(metfile,&CARDADATA);
*/




FREE_DATA_STRUCT(DATA);





return 0;

}
