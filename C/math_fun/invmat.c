#include<stdio.h>
#include <malloc.h>
int invmat(double **MAT,double **INVMAT,int n)
{
    
    double ratio,a;
    int i, j, k;


printf("MAT = %i\n",MAT);
printf("MAT[0] = %i\n",MAT[0]);
printf("MAT[0][0] = %f\n",MAT[0][0]);


printf("IMAT = %i\n",INVMAT);
printf("IMAT[0] = %i\n",INVMAT[0]);
printf("IMAT[0][0] = %f\n",INVMAT[0][0]);


printf("OK done - now working it out\n");


    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i==(j))
                INVMAT[i][j] = 1.0;
            else
                INVMAT[i][j] = 0.0;
        }
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i!=j){
                ratio = MAT[j][i]/MAT[i][i];
                for(k = 0; k < n; k++){
                    MAT[j][k] -= ratio * MAT[i][k];
                    INVMAT[j][k] -= ratio * INVMAT[i][k];
                }
            }
        }
    }
    for(i = 0; i < n; i++){
        a = MAT[i][i];
        for(j = 0; j < n; j++){
            MAT[i][j] /= a;
            INVMAT[i][j] /= a;
        }
    }
    return 0;
}
