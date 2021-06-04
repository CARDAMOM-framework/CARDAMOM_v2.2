#pragma once
double meanptr(double *P,int ne){
int n;
double meanP;
for (n=0;n<ne;n++){meanP=meanP+P[n];}
return meanP/(double)ne;



}





