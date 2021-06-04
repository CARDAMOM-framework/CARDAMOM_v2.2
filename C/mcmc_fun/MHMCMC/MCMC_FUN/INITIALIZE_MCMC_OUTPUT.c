#pragma once
int INITIALIZE_MCMC_OUTPUT(PARAMETER_INFO PI,MCMC_OUTPUT *MCOUT, MCMC_OPTIONS MCO){

/*defining MCMC output structure*/
MCOUT->best_pars=calloc(PI.npars*MCO.nchains,sizeof(double));
/*set to zero, will change to 1 when MCMC is complete!*/
MCOUT->complete=0;


return 0;


}

