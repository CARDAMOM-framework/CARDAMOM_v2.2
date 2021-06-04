typedef struct{
/*Adapt step size every N iterations*/
int nADAPT;
/*Adapt step size for XX fraction of full run*/
double fADAPT;
/*Number of requested parameter vectors*/
int nOUT;
/*Print info every N solutions (set to 0 for silent)*/
int nPRINT;
/*Write every N solutions*/
int nWRITE;
/*Append previous file with same options (1-0)*/
int APPEND;
/*output file name*/
char outfile[200];
/*saved step filename*/
char stepfile[200];
/*saved step filename*/
char startfile[200];
/*return best fit*/
int returnpars;
/*random ini pars*/
int randparini;
/*int fixedpars (1-0)= fix initial values that are not equal to -9999*/
int fixedpars;
/*minstepsize*/
double minstepsize;
/*number of chains*/
int nchains;
/*mcmc type*/
int mcmcid;
}MCMC_OPTIONS;


/*defining all immediately useful MCMC OUTPUTS
 * such as best fitting parameters
 * etc.*/
typedef struct{
/*1 or 0*/
int complete;
/*add any stats you like*/

/*best parameters here*/
double *best_pars;
}MCMC_OUTPUT;


/*defining all counters*/
typedef struct{
/*number of accepted solutions*/
int ACC;
/*number of recently accepted solutions*/
int ACCLOC;
/*local acceptance rate*/
double ACCRATE;
/*number of iterations*/
int ITER;
/*running normalized parameter mean*/
double *parmean;
/*running standard deviation*/
double *parstdev;
/*running covariance*/
double **parcovariance;
/*Cholesky decomp mat*/
double **parcholesky;
/*stepsize amplitude*/
double amp;
/*running covariance sample count*/
int covsamplecount;
}COUNTERS;



/*PARAMETER STRUCTURE DEFINED HERE*/
/*here all parameter information is defined*/


typedef struct{
/*number of parameters*/
int npars;
/*maximum parameter values*/
double *parmax;
/*minimum parameter values*/
double *parmin;
/*initial parameter values*/
double *parini;
/*parfix*/
double *parfix;
/*initial step size*/
double *stepsize;
/*PCArotation*/
double **PCArotation;
/*circular/log/etc parameters*/
int *transform;
}PARAMETER_INFO;

