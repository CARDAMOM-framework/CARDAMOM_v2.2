#pragma once
#include <math.h>
/*CONVERTING PARAMETERS TO 0-1 range*/


double par2nor(double p, double mn, double mx){
return log(p/mn)/log(mx/mn);}


/*AND VISE VERSA*/


double nor2par(double p, double mn, double mx){
return mn*pow(mx/mn,p);}
