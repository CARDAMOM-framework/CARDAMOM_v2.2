
#pragma once
/*This function is called by DALEC_FIRE.c DALEC_CDEA.c DALEC_BUCKET.c, etc.*/

double offset(double const L, double const w)
{

/*solution to t = f*exp(-t^2)*/
/*see dalecstepfunction2*/

double mxc[7]={0.000023599784710, 0.000332730053021,    0.000901865258885,  -0.005437736864888,  -0.020836027517787,   0.126972018064287,   -0.188459767342504};

double lf=log(L-1);
double os=mxc[0]*pow(lf,6) + mxc[1]*pow(lf,5) + mxc[2]*pow(lf,4) + mxc[3]*pow(lf,3) + mxc[4]*pow(lf,2) + mxc[5]*lf +mxc[6];

os=os*w;

return os;
}


