#pragma once
int WRITE_RESULTS(double *PARS, double PROB,PARAMETER_INFO PI,MCMC_OPTIONS MCO){

int n;



FILE *fileout=fopen(MCO.outfile,"ab");
/*FILE *filestep=fopen(MCO.stepfile,"wb");*/
for (n=0;n<PI.npars;n++){
        fwrite(&PARS[n],1,sizeof(double),fileout);
       	/*fwrite(&PI.stepsize[n],1,sizeof(double),filestep);*/
}
    
/*writing likelyhood*/
/*NOTE: As of July 11th 2014, probability no longer written to file*/
/*Probability is a "re-derivable" quantity, therefore if needed, it can
either (a) be re-derived using the MODEL_LIKELIHOOD function, or (b) 
written to a separate file.*/
/*fwrite(&PROB,1,sizeof(double),fileout);*/


fclose(fileout);


return 0;


}
