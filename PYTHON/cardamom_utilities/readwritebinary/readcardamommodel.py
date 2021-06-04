#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:41:18 2020

@author: gquetin

Read data in the CARDAMOM model structure

Take function from analycbr

"""

import numpy as np
import os
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 


import library_cbr as lcbr
imp.reload(lcbr)

def cbflabels():
    
    CBFmetlabels = ('Time [Days since Jan 01, 2001]','min temp [C]', 'max temp [C]',
                'Solar irradiance [MJ/m2/day]','CO2','Day of year',
                'Burned area [m2/m2]', 'VPD [hPa]','Precip. [mm/day]')
    
    return CBFmetlabels


def library_parameter(modelid,modeldir):
    parameterfile = 'projects/CARDAMOM_MODELS/DALEC/' + 'DALEC_' + str(modelid) + '/'+'PARS_INFO_'+str(modelid)+'.c'
    A = open(modeldir + parameterfile).readlines()
    
    parnames = {}
    for idx,a in enumerate(A):
        if 'CARDADATA' in a and 'CARDADATA' in A[idx+1]:
            
            if not np.any([a[0:2] == '//',A[idx+1][0:2] == '//']):
                b1 = a.find('[')
                b2 = a.find(']')
                pmin1 = a.find('=')
                pmin2 = a.find(';')
                pmax1 = A[idx+1].find('=')
                pmax2 = A[idx+1].find(';')
                parnum_here = int(a[b1+1:b2])
                
                bb = A[idx-(3-np.max([idxx for idxx, aa in enumerate(A[idx-4:idx+1]) if aa.strip() =='']))]
                n1 = bb.find('/*')
                n2 = bb.find('*/')
                parnames[parnum_here] = {'name':bb[n1+2:n2],
                                         'parmin':a[pmin1+1:pmin2],
                                         'parmax':A[idx+1][pmax1+1:pmax2]}

           
    return parnames





def get_parnames(modelid,output='shortnames',modeldirshort='/repos/models/cardamom/C/'  ):
    modeldir = MACHDIR+modeldirshort
    param_dict_org = library_parameter(modelid,modeldir)
    pnames_org = [param_dict_org[pn]['name'] for pn in param_dict_org]
    
    
    
    
    
    if output == 'shortnames':
        cbr_parsname_dict = lcbr.cbr_par_dictionary()
        pnames_org_short = [cbr_parsname_dict[pn][0] for pn in pnames_org]
        return pnames_org_short
    elif output == 'longnames':
        return pnames_org
    elif output == 'dictionary':
        return param_dict_org
    else:
        print('Not Available Yet')
        return