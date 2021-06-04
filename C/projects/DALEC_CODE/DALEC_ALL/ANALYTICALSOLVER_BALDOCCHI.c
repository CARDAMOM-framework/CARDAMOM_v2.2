#pragma once

double ANALYTICALSOLVER_BALDOCCHI(double const *met_list, double const *var_list, double const *wc_adeb_list, double const *wj_adeb_list, double const Rd, double const RH, double const ROOT_CALC)
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
    
    
    /* wj_adeb_list[0] = aj */
	/* wj_adeb_list[1] = dj */
	/* wj_adeb_list[2] = ej */
	/* wj_adeb_list[3] = bj */
	/* wc_adeb_list[0] = ac */
	/* wc_adeb_list[1] = dc */
	/* wc_adeb_list[2] = ec */
	/* wc_adeb_list[3] = bc */
	
    /******************************************/
	double aj, dj, ej, bj, ac, dc, ec, bc;
	
	aj = (double)wj_adeb_list[0];
	dj = (double)wj_adeb_list[1];
	ej = (double)wj_adeb_list[2];
	bj = (double)wj_adeb_list[3];
	ac = (double)wc_adeb_list[0];
	dc = (double)wc_adeb_list[1];
	ec = (double)wc_adeb_list[2];
	bc = (double)wc_adeb_list[3];


    /***********************************************/
    /*ANALYTICAL SETUP*/

    /*initialize analytical setup variables*/
    double alpha, beta, sgamma, theta1;
    double pj, qj, rj, pc, qc, rc, Assim_j, Assim_c, GPP;
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
    
    /* Root 1*/
    if (ROOT_CALC == 1){
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 4.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 4.*pi) /3. ) - pc/3.;
	
    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
    }
    
    
    if (ROOT_CALC == 2){
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 2.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 2.*pi) /3. ) - pc/3.;
    
    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
    }
    
    
    if (ROOT_CALC == 3){
    Assim_j=(double)-2.*pow((pow(pj,2.)-3.*qj)/9.,0.5) * cos( (acos((2.*pow(pj,3.) - 9.*pj*qj + 27.*rj)/(2.*pow((pow(pj,2.) - 3.*qj),1.5))) + 0.*pi) /3. ) - pj/3.;
    Assim_c=(double)-2.*pow((pow(pc,2.)-3.*qc)/9.,0.5) * cos( (acos((2.*pow(pc,3.) - 9.*pc*qc + 27.*rc)/(2.*pow((pow(pc,2.) - 3.*qc),1.5))) + 0.*pi) /3. ) - pc/3.;
    
    /* pick the smaller A to be GPP */
    /* GPP UNITS ARE umol m^-2 s^-1 */
    GPP = (Assim_j<=Assim_c) ? Assim_j : Assim_c;
    }
    
    
    return GPP;
    
    }