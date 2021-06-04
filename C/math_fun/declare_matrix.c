#pragma once
#include <stdlib.h>

int declare_matrix( double ***M, int num_rows, int num_columns)
{
    
    double *pool;
    double *curPtr;
    //(step 1) allocate memory for array of elements of column

    M[0] = (double** )calloc(num_rows, sizeof(double* ));

    //(step 2) allocate memory for array of elements of each row
    pool = (double *)calloc( num_columns * num_rows, sizeof(double));

    // Now point the pointers in the right place
    curPtr = pool;
int i;
    for(i = 0; i < num_rows; i++)
    {
        *(M[0] + i) = curPtr;
         curPtr += num_columns;
    }

return 0;

}
