
#include "../../auxi_fun/okcheck.c"
#include "math.h"
#include <stdio.h>
#include "CARDAMOM_READ_BINARY_DATA.c"
/*syntax CARDAMOM_READ_BINARY_CARDADATA(char *filename,CARDADATA *CARDADATA)*/




/* the gateway function */
/*NOTE: this function can be used with any valid model ID (stored in CARDADATA.ID)*/
int main(int argc, char * arg[])
{
 
/*declaring loop variable n*/ 

/*storing command line inputs as 2 files*/



/*declaring data structure*/
DATA DATA;

/*Runs all model info functions to assemble models and generate outputs*/ 
int id=0,cml;
for (id=0;id<10000;id++){
DATA.ID=id;
DATA.assemble_model=0;

CARDAMOM_MODEL_LIBRARY(&DATA);

if (DATA.assemble_model==1){WRITE_MODEL_STACKS_TO_FILE(&DATA,arg[1]);}
}

/*make file with output summary*/



/*
CARDAMOM_READ_BINARY_DATA(metfile,&CARDADATA);
*/




FREE_DATA_STRUCT(DATA);





return 0;

}
