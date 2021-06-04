
#pragma once
#include "../DALEC_CODE/DALEC_ALL/DALEC_MODULE.c"

int WRITE_MODEL_STACKS_TO_FILE(DATA * DATA,char * pathname){
char fname[100], fullfile[500];
strcpy(fullfile,pathname);
sprintf(fname,"/projects/CARDAMOM_MODELS/DALEC/DALEC_%i/dalec_%i_pars.txt",DATA->ID, DATA->ID);
strcat(fullfile,fname);
printf("****\n");
printf("file = %s\n",fullfile);

FILE *fdf;
fdf=fopen(fullfile,"w");
filediag(fdf,fullfile);

int n;

DALEC *MODEL=(DALEC *)DATA->MODEL;

fprintf(fdf, "DALECmodel.nopars=%i;\n",MODEL->nopars);
fprintf(fdf, "DALECmodel.nopools=%i;\n",MODEL->nopools);
fprintf(fdf, "DALECmodel.nomet=%i;\n",MODEL->nomet);
fprintf(fdf, "DALECmodel.nofluxes=%i;\n",MODEL->nofluxes);



for (n=0;n<DATA->nopars;n++){
fprintf(fdf,"DALECmodel.parname{%i}='%s';\n",n+1,DATA->parname[n]);
fprintf(fdf,"DALECmodel.parmin(%i)=%8.4f;\n",n+1,DATA->parmin[n]);
fprintf(fdf,"DALECmodel.parmax(%i)=%8.4f;\n",n+1,DATA->parmax[n]);}

fclose(fdf);

printf("writing success!!!");
return 0;}
