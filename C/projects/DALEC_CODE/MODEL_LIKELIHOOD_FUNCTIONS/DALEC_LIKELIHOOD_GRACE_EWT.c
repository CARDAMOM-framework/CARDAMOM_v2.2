double DALEC_LIKELIHOOD_GRACE_EWT(DATA D){
/*Data structure, includes model and data*/

/*EWT constraint*/
double tot_exp=0;
double mewt=0, mewtm=0;
int n,dn;
double P=0;
if (D.newt>0){
/*Relative volumetric constraint*/

for (n=0;n<D.newt;n++){

/*Note: constraint imposed on mean between t and t+1 of H2O pools*/
dn=D.ewtpts[n];

/*mean model EWT*/
mewtm+=(D.M_POOLS[D.nopools*dn+6] + D.M_POOLS[D.nopools*dn+7]+
D.M_POOLS[D.nopools*(dn+1)+6] + D.M_POOLS[D.nopools*(dn+1)+7]
)/(D.newt*2);
/*Mean obs ewt - already aligned between t and t+1*/
mewt+=(D.EWT[dn])/(D.newt);

}

for (n=0;n<D.newt;n++){dn=D.ewtpts[n];tot_exp+=pow((
(D.M_POOLS[D.nopools*dn+6]+ D.M_POOLS[D.nopools*dn+7]
+D.M_POOLS[D.nopools*(dn+1)+6]+ D.M_POOLS[D.nopools*(dn+1)+7])/2
-mewtm - D.EWT[dn]+mewt)/D.ewt_obs_unc,2);}
P=P-0.5*tot_exp;}



/*Returning log-probability*/
return P;
}
