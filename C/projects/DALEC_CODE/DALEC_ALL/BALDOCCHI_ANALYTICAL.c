#pragma once
#include "./ANALYTICALSOLVER_BALDOCCHI.c"

double* BALDOCCHI_ANALYTICAL(double const *met_list, double const *var_list)
{
    /****************** INPUTS *****************/
    /* met_list[0] = mintemp */
    /* met_list[1] = maxtemp */
    /* met_list[2] = co2 */
    /* met_list[3] = rad */
    /* met_list[4] = VPD */

    /* var_list[0] = m */
    /* var_list[1] = b */
    /* var_list[2] = gb */
    /* var_list[3] = Jmax */
    /* var_list[4] = Vmax */
    /* var_list[5] = LAI */
    /* var_list[6] = k */
    /* var_list[7] = c1 */
    /* var_list[8] = c2 */
    /******************************************/

    /*initialize intermediate variables*/
    double Oi, psfc, T, RH, tau, gamma, Kc, Ko;
    double cold_inhibition, heat_inhibition, Vm;
    double APAR, a, b, c, J1, J2, J, Rd;
    double GPP, GPP_root1, GPP_root2, GPP_root3;
    double aj, dj, ej, bj, ac, dc, ec, bc;
    double Wm2_umol;

	Wm2_umol = (double)2.02; /*conversion of solar radiation to photosynthetically active radiation, accounts for only ~45% light in visible, Mavi & Tupper (2004) */
    psfc=(double)101325.; /* [Pa]*/
    Oi=(double)0.209*psfc; /*oxygen*/
    T=(double)0.5*(met_list[1]+met_list[0]);
    RH=(double)1.-met_list[4]/(6.1094*exp( (17.625*T) / (T+243.04) )); /*relative humidity, <1*/

    tau=(double)2600.*pow(0.57,((T-25.)/10.));
    gamma=(double)Oi/(2.*tau); /*co2 compensation point*/

    Kc=(double)30.*pow(2.1,((T-25.)/10.)); /*michaelis constant for carboxylation*/
    Ko=(double)30000.*pow(1.2,((T-25.)/10.)); /*michaelis constant for oxidation*/

    /*compute temperature-adjusted carboxylation rate*/
    cold_inhibition=(double)1+exp(0.25*(10.-T));
    heat_inhibition=(double)1+exp(0.4*(T-40.));
    Vm=(double)var_list[4]*pow(2.1,((T-25.)/10.))/(cold_inhibition*heat_inhibition);

    /*solve quadratic to get J value*/
    APAR=(double)met_list[3]*(1-exp(-var_list[6]*var_list[5]))*1e6/86400.*Wm2_umol; /*trailing constants here ensure translation from MJ m^-2 d^-1 to umol m^-2 s^-1 PAR, which is required for our analytical setup*/
    a=(double)0.7;
    b=(double)-(var_list[3]+0.385*APAR);
    c=(double)0.385*var_list[3]*APAR;
    J1=(double)(-b+pow((pow(b,2.)-4*a*c),0.5)) / (2*a);
    J2=(double)(-b-pow((pow(b,2.)-4*a*c),0.5)) / (2*a);
    J = (J1<=J2) ? J1 : J2;

    /*compute temperature-adjusted leaf respiration*/
    Rd=(double)0.015*Vm*pow(2.4,((T-25.)/10.)) / (1.+exp((1.3*(T-55.))));

    /*set the a,d,e,b parameters for wj*/
    aj=(double)J;
    dj=(double)gamma;
    ej=(double)4.5;
    bj=(double)10.5*gamma;

    /*set the a,d,e,b parameters for wc*/
    ac=(double)var_list[4];
    dc=(double)gamma;
    ec=(double)1.;
    bc=(double)Kc*(1.+Oi/Ko);

	/**************************/
	/* check dGPP/dCa */
	double wj_adeb_list[] = {aj,dj,ej,bj};
	double wc_adeb_list[] = {ac,dc,ec,bc};
	
	double ROOT_CALC, GPP_Ca_high, GPP_Ca_low, dGPPdCa;
	double met_list_high[] = {met_list[0],met_list[1],40.,met_list[3],met_list[4]}; /* Switch in high value atmospheric CO2 */
	double met_list_low[] = {met_list[0],met_list[1],30.,met_list[3],met_list[4]}; /* Switch in high value atmospheric CO2 */
	ROOT_CALC = 1.;
	
	
	GPP_Ca_high = ANALYTICALSOLVER_BALDOCCHI(met_list_high, var_list, wc_adeb_list, wj_adeb_list, Rd, RH, ROOT_CALC);
	GPP_Ca_low = ANALYTICALSOLVER_BALDOCCHI(met_list_low, var_list, wc_adeb_list, wj_adeb_list, Rd, RH, ROOT_CALC);
	dGPPdCa = GPP_Ca_high - GPP_Ca_low;
	


    /***********************************************/
    /*ANALYTICAL SETUP*/

    /*initialize analytical setup variables*/
    double alpha, beta, sgamma, theta1;
    double pj, qj, rj, pc, qc, rc, Assim_j, Assim_c;
    double pi=3.1415927;

    /*compute combined variables*/
    alpha=(double)1.647 + var_list[1]/var_list[2] - var_list[0]*RH;
    beta=(double)met_list[2]*(var_list[2]*(var_list[0]*RH-1.647) - 2.*var_list[1]);
    sgamma=(double)pow(met_list[2],2.)*var_list[1]*var_list[2];
    theta1=(double)var_list[2]*var_list[0]*RH - var_list[1];

    /*compute coefficients for cubic (wj)*/
    pj=(double)(ej*beta + bj*theta1 - aj*alpha + ej*alpha*Rd) * (1/(ej*alpha));
    qj=(double)(ej*sgamma + (bj*sgamma)/met_list[2] - aj*beta + aj*dj*theta1 + ej*Rd*beta + Rd*bj*theta1) * (1./(ej*alpha));
    rj=(double)(-aj*sgamma + (aj*dj*sgamma)/met_list[2] + ej*Rd*sgamma + (Rd*bj*sgamma)/met_list[2]) * (1./(ej*alpha));

    /*compute coefficients for cubic (wc)*/
    pc=(double)(ec*beta + bc*theta1 - ac*alpha + ec*alpha*Rd) * (1/(ec*alpha));
    qc=(double)(ec*sgamma + (bc*sgamma)/met_list[2] - ac*beta + ac*dc*theta1 + ec*Rd*beta + Rd*bc*theta1) * (1./(ec*alpha));
    rc=(double)(-ac*sgamma + (ac*dc*sgamma)/met_list[2] + ec*Rd*sgamma + (Rd*bc*sgamma)/met_list[2]) * (1./(ec*alpha));

    /*compute assimilation rates for wj and wc separately, according to analytical solution*/
    /* Update 04.07.2020 by G. Quetin to calculate multiple roots and switch */
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 4.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 4.*pi) /3. ) - pc/3.;

    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP_root1 = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
	
	if (dGPPdCa > 0.){
	GPP = GPP_root1;
	/*printf("Root1: %f\n",GPP);*/
	}
	
	/* Calculate the second root of the polynomial fit */
	if ((dGPPdCa <= 0.) & (GPP_root1 < 0)){
	/*compute assimilation rates for wj and wc separately, according to analytical solution*/
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 2.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 2.*pi) /3. ) - pc/3.;

    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP_root2 = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
    GPP = GPP_root2;
    /*printf("Root2: %f\n",GPP);*/
	}
	
	/* Calculate the third root of the polynomial fit */
	if ((dGPPdCa <= 0.) & (GPP_root1 > 0)){
	/*compute assimilation rates for wj and wc separately, according to analytical solution*/
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 2.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 2.*pi) /3. ) - pc/3.;

    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP_root3 = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
    GPP = GPP_root3;
    /*printf("Root3: %f\n",GPP);*/
	}
	
	/* Do not allow negative values of GPP */
	if (GPP <= 0){
	GPP = 0;
	}
	
	
	/*printf("Good after GPP leaf level");*/

    /***********************************************/
    /*SCALE TO CANOPY LEVEL*/

    double GPP_scaled; /*E0*/
    /*E0=(double)(var_list[7]*pow(var_list[5],2.))/(var_list[8] + pow(var_list[5],2.));
    GPP_scaled=(double)(E0*met_list[3]*1e6/0.327/86400.*GPP)/(E0*met_list[3]*1e6/0.327/86400. + GPP);*/
    GPP_scaled=(double)(GPP* (1-exp(-var_list[6]*var_list[5]))/var_list[6] );

    /***********************************************/
    /*CALCULATE WATER USE EFFICIENCY*/

    double Cs, gs, Gs_canopy;//, iWUE;
    Cs=(double)met_list[2]*10.-GPP/var_list[2]; /*Cs is in umol mol^-1*/
    gs=(double)var_list[0]*GPP*RH/Cs+var_list[1]; /*gs is in mol m^-2 s^-1*/
    Gs_canopy=(double)(gs* (1-exp(-var_list[6]*var_list[5]))/var_list[6] );
    //iWUE=(double)GPP_scaled/Gs_canopy; /*iWUE is in umol mol^-1*/
    //iWUE=(double)iWUE/10/100.*12.0107*1000./18.01528; /*iWUE now in gC/kgH2o *hPa, same as the original IWUE parameter (optimized)*/

    /***********************************************/
    /*RETURN AN ARRAY*/

    GPP_scaled = GPP_scaled*1.03775; /*convert from umol m^-2 s^-1 to gC m^-2 d^-1*/

    //iWUE=(double)GPP_scaled/(Gs_canopy*18.01528*86400./1000.); /*iWUE is in gC / kgH2O

    static double return_arr[2]; /*need static variable to be able to return array?*/
    return_arr[0]=GPP_scaled;
    //return_arr[1]=iWUE;
    return_arr[1]=Gs_canopy*1.6*18.01528*86400./1000.; //units of kg H2O m^-2 day^-1

  return return_arr;
}
