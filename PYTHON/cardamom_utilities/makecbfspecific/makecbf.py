#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:13:21 2018

@author: gquetin

Make a .cbf file for input to CARDAMOM.c

Inputs:
    1. Other cbf files (for met, fire, CO2 etc., binar or .mat)
    2. Other sources of fire (climatologies, LENS, empirical)
    3. CO2 (experiment sources etc.)
    4. New met (CMIP5 netcdfs, LENS netcdfs, reanalysis) [consider a preprocessing]
    
Outputs:
    1. cbf (binary)
    2. cbf (netcdf, optional)
    3. cbf (.mat, optional)
    4. text log (.txt, optional)
    
Options to make a set of dates and data sources?

"""


import numpy as np

import os
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 

import cardamom_map_utility as cardmap
import readwritebinary as rwbin

imp.reload(rwbin)
imp.reload(cardmap)


def cbflabels(): # putting in readcardamommodel.py
    
    CBFmetlabels = ('Time [Days since Jan 01, 2001]','min temp [C]', 'max temp [C]',
                'Solar irradiance [MJ/m2/day]','CO2','Day of year',
                'Burned area [m2/m2]', 'VPD [hPa]','Precip. [mm/day]')
    
    return CBFmetlabels

CBFmetlabels = cbflabels()

def selectcbf(rowcol,cbffilelist):
    '''
    Take a row and column input (string), find and open the matching
    cbf file in the directory.
    '''
            
    lat, lon = cardmap.rowcol_to_latlon([rowcol])[0]
    dslat, dslon = cardmap.lon180_to_lon360(lat,lon)

    cbffilename = None
    for cbff in cbffilelist:
        if rowcol in cbff:
            cbffilename = cbff
            break
    
    if cbffilename == None:
        print('ERROR: No matching cbf file for point ' + rowcol)
    
    CBF = rwbin.read_cbf_file(cbffilename)
    
    return CBF




def datetime_from_cbf(dayarray,unitt = 'D',refdate='2001-01-01'):

    dayssince = np.array([np.timedelta64(int(s), unitt) for s in dayarray])
    
    cbftime = np.datetime64(refdate,unitt) - np.timedelta64(int(1), unitt) + dayssince
    cbftimens = np.array([np.datetime64(ts,'ns') for ts in cbftime])
    
    return cbftimens




def cardamom_obsnames(version = 'v1'):
    
    if version == 'v1':
        obsnames = ['sif','lai','nbe']
    else:
        print('no such version: ' + version)
        obsnames = []
    
    return obsnames



def cardamom_othernames():
    
    #    CBF['OBS']:
#        The three following columns are required
#        Column 1 = GPP (or SIF, in case of SIF make sure CBF.gppabs=0)
#        Column 2 = LAI (mean LAI is constrained 
#        Column 3 = NEE (annual NEE is constrained if CBF.neeiav=1)
#    
#        The three following columns are optional
#        Column 4 = Biomass
#        Column 5 = ET 
#        Column 6 = GRACE equivalent water thickness
#        
#    CBF['OTHERPRIORS']:
#        Numbers per python indexing
#        0: Biomass Constraint
#        1:
#        2: Mean Fire Emission [gC/m2/day]
#        3:
#        4: Mean LAI constraint
#            if individual LAI values are prescribed at some time steps in CBF.OBS, 
#            CBF.OTHERPRIORSUNC(5) overrides these but uses the time steps to 
#            calculate mean modeled LAI. If no LAI values are prescribed in 
#            CBF.OBS, then meanmodeled LAI is calculated based on all model time steps.
#        5: Mean GPP [gC/m2/day]
#        
#    CBF['OTHERPIRORSUNC']:
#        0: Biomass Uncertainty Factor
#        1:
#        2: Mean fire emssions (for absolute uncertainty prescription, prescribe as -ve value)
#        3:
#        4: Mean LAI uncertainty factor
#        5: Mean GPP uncertainty factor (for absolute uncertainty prescription, prescribe as -ve value)
#        
#    CBF['PARPRIORS']:
#        
#    CBF['PARPRIORUNC']:
# obsnames={'GPP','LAI','NBE','ABGB','ET','EWT','BAND1','BAND2','BAND3','BAND4','SOM'};
    
    CBFlegend = {}
    CBFlegend['OBS'] = ['gpp',     # (or SIF, in case of SIF make sure CBF.gppabs=0)
                        'lai',     # (mean LAI is constrained )
                        'nbe',     # (annual NBE is constrained if CBF.neeiav=1)
                        'tbiomass', # Optional
                        'et',      # Optional
                        'grace',
                        'band1',
                        'band2',
                        'band3',
                        'band4',
                        'som']   # Optional, equivalent water thickness

        
    CBFlegend['OTHERPRIORS'] = ['biomass',
                                  'None',
                                  'fire', # Mean Fire Emission [gC/m2/day]
                                  'None',
                                  'meanlai', # Mean LAI constraint
        #                                      if individual LAI values are prescribed at some time steps in CBF.OBS, 
        #                                      CBF.OTHERPRIORSUNC(5) overrides these but uses the time steps to 
        #                                      calculate mean modeled LAI. If no LAI values are prescribed in 
        #                                      CBF.OBS, then meanmodeled LAI is calculated based on all model time steps.
                                  'meangpp'  # [gC/m2/day]
                                  ] +['None']*44

        
    CBFlegend['OTHERPRIORSUNC'] = ['biomass_unc',
                                      'tbiomass_unc', # Uncertainty for biomass time series
                                      'fire_unc', # Mean Fire Emission [gC/m2/day]
                                      'None',
                                      'meanlai_unc', # Mean LAI constraint
                                      'meangppz_unc',# [gC/m2/day]
                                      'None',
                                      'som_unc',
                                      ] +['None']*42

        
    CBFlegend['PARPRIORS'] = ['par'+str(jj) for jj in range(0,50)]
        
    CBFlegend['PARPRIORUNC'] = ['parunc'+str(jj) for jj in range(0,50)]
    
    return CBFlegend 


# %%

if __name__ == '__main__':
    
    print('Main file functionality has been moved to run-makecbf.py')
            