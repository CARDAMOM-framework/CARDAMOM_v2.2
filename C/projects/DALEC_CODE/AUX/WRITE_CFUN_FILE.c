void WRITE_CFUN_FILE(char *site,double *CFUN, int N)

{

char filename[]="CFUN_";
strcat(filename,site);

FILE *file;

file=fopen(filename, "a+");

int n;
for (n=0;n<N;n++){
fprintf(file,"%f \n",CFUN[n]);}


fclose(file);




}
