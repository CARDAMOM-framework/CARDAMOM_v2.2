#pragma once
/*This function contains to model-specific structure*/
/*starting with constant terms showing up across scripts*/
/*unclear where this should be declared, but "pragma once" solves problem for now*/

/*EDCdiagnostic structure*/
/*While this structure type is declared as part of the "DALEC" structure, should include here*/
struct EDCDIAGNOSTIC{
int EDC;
double pEDC;
int DIAG;
/*allocating space for 100 checks (in case future additions are made, etc.).*/
int PASSFAIL[100];
/*SWITCH determines which EDCs are tested*/
/*SWITCH values are read from DATA.otherpriors*/
/*SWITCH values that are not 1 (e.g. -9999) are assumed to be zero*/
int SWITCH[100];
/*EDCPROB: this term is a 0 - 1 assessment (currently qualitative) to determine
the viability of individual EDCs*/
double EDCPROB[100];
/*nedc is the number of EDCs to be tested*/
int nedc;
/*Temporary value - EQF - steady state proximity factor*/
/*Currently stored in OTHERPRIORUNC (39:42)*/
double EQF;};


typedef struct DALEC{
int nopools;
int nopars;
int nofluxes;
int nomet;
int (*dalec)(DATA,const double *);
int (*edc1)(const double *, DATA, struct EDCDIAGNOSTIC * EDCD);
int (*edc2)(const double *, DATA, struct EDCDIAGNOSTIC * EDCD);
/*contains all the EDCD relevant info*/
struct EDCDIAGNOSTIC * EDCD;
}DALEC;
