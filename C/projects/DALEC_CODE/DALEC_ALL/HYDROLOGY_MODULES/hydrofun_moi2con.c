double hydrofun_moi2con(double sm,double k0, double b){


//m3/m3 to m/s
double moi2con=pow(k0*sm,2*b+3);


return moi2con;
}
