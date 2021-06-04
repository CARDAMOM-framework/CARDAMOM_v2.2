#pragma once


/*mean matrix from double pointer routine*/
double mean_pool(double *PA,int p,int nc, int nopools){
/*declarations*/
int r,c;
double meanpool=0;
/*deriving mean of pool p*/
for (c=0;c<nc;c++){
meanpool=meanpool+PA[c*nopools+p]/(double)nc;
}
/*returing meanpool value*/
return meanpool;}






