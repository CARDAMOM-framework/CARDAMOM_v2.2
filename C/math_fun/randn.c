#pragma once
double randn(){

double pi=3.141592653589793;
double r1=(double)random()/RAND_MAX;
double r2=(double)random()/RAND_MAX;


double rn=sqrt(-2*log(r1)) * cos(2*pi*r2);

return rn;

}

