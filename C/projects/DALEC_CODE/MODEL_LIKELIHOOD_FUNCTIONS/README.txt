******UPDATES*******

-DALEC_MLF.c will become the standard MLF function.
-Other MLF.c functions will eventually be decommissioned
-Decommissioning for any MLF function possible by model declaration for compatibility with DALEC_MLF.c
-See ID=804 and later models for example


******GENERAL COMMENTS*******

These functions are to be used for data assimilation projects.
As these are generally similar these are all kept in the same folder.
Any functions used in two or more models are kept in DALEC_ALL_LIKELIHOOD.c
For the time being, similar functions for different models (e.g. EDC functions, MLF functions) are syncronized (as it is simpler than making a generic function at this stage). 
HOWEVER: any subroutine that is called by two or more functions is kept in a common space.

*****SPECIFICS************

The first position (0) of DATA.otherpriors is used to define total live biomass.
The 31st-50th positions (30-49) of DATA.otherpriors are used to define EDC switches (only applicable if EDC diagnostic mode is set to 2).
