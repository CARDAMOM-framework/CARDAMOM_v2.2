#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:10:27 2019

@author: gquetin

Convert a .cbr file to a netcdf file, or a set of .cbr files to one netcdf with
chain information

"""

import scipy.io as sio
import numpy as np
import xarray as xr

import os
import glob
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 
import file_organization_utility as forguti

import readwritebinary as rwbin
import readcardamommodel

imp.reload(forguti)
imp.reload(rwbin)



def cbr2netcdf(parfilelist,cutfrac,varnames,INFO):
    '''
    Convert a .cbr file or list of .cbr files into a netcdf with parameter names
    '''
    Nlist = []
    pars = []
    for parfile in parfilelist:
        if os.path.exists(parfile): # and os.path.getsize(parfile) == filesize:
            parsraw = rwbin.read_cbr_file(parfile,INFO)
            N = int(cutfrac*parsraw.shape[0])
            pars.append(parsraw[N:,:])
            Nlist.append(parsraw.shape[0] - N)
    
    
    pars = np.concatenate(pars,axis=0)   
    
   

    # Vector of parameters sets, just labeled 1 to 1000
    diffpars = np.arange(1,pars.shape[0]+1,1) #number of models run with differen parameters
    


    #
    # Create dataset of outputs
    dspar = xr.Dataset()
    
    for idx, pl in enumerate(varnames):
        dspar[pl] = xr.DataArray(pars[:,idx],coords=[diffpars],
               dims = ['parameterset'])
#                    dspar[pl].attrs = {'longname':poolnames(output='dict')[pl],
#                                       'units':poolnames(output='dictunits')[pl]}
        
    dspar.attrs['chainnums'] = np.cumsum([0]+Nlist)
    
    return dspar




# %%

if __name__ == '__main__':
    
    plotsave = True
    datasave = True
    mapplot = False
    makeplots = False
    
    # Set up folders based on the home folder
    datadir, dataoutdir = forguti.folderdeff(MACHDIR)
    
    exper = 'p16_long' #'p16_long' #'p6_uncmatch_v1_smthnbe'
    model = 'cruncep_osse5' #'erainterim'
    cutfrac = .75 # may have this labeled wrong, says how much to keep rather than throw away
    filesize = 144000 # 144000
    
    dircbr = MACHDIR + '/Google Drive/DATA/DALEC/'+exper+'/'+model+'/cbr/'
    plotdir = datadir + 'Analysis/DALEC_CMIP5/p16_long/'+'OSSE1plus' +'/'
    
   
    modelid = 821 
    
    
    pnames_org_short = readcardamommodel.get_parnames(modelid,output='shortnames')
            
    INFO = {'nopars':len(pnames_org_short),
            'latterhalf':0}

    lparnames = pnames_org_short 
    
    import cardamom_output2netcdf as carout2nc
    filelist = glob.glob(dircbr + '*.cbr')
    
    wildcard = list(set([carout2nc.file_to_colrow(fl)[1] for fl in filelist]))
    ptlist = list(set([carout2nc.file_to_colrow(fl)[0] for fl in filelist]))
    
    itrpick = '100M'
    dict_chain = {'100M':[1,2,3,4,5,6,7,8],
                  '20M':[9,10,11,12,13,14,15,16]}
    
    chainlist = dict_chain[itrpick]
    
    #%%
    
    
    dict_par = {}
    for pt in ptlist[0:]:
        
        
        parfilelist = []
        for ch in chainlist[0:]:
        
            parfilelist.append('{}{}_{}_{}.cbr'.format(dircbr,wildcard[0],pt,ch))
            

                
        dspar = cbr2netcdf(parfilelist,cutfrac,varnames = pnames_org_short,INFO=INFO)
        dict_par[pt] = dspar
                
    #%%
    
    plotcheck = True
    if plotcheck == True:
        
        
        import matplotlib.pyplot as plt
        
        for pt in ptlist[0:]:
            
            dsplt = dict_par[pt]
            pltlist = list(dsplt.data_vars)
            fig, axarr = plt.subplots(int(len(pltlist)/6),6,figsize=(20,8))
            for jjv, var in enumerate(pltlist):
                ax = axarr.flatten()[jjv]
                ax.hist(dsplt[var],bins=40)
                ax.set_title(pt + ' : ' + var)
            plt.tight_layout()
            
            if plotsave is True:    
                plotname = 'cbrpt_parhist'
                                
                savetag = "_".join([pt,model,itrpick])
                fig.savefig(plotdir + plotname + savetag)
                #plt.savefig(plotdir + plotname + savetag +'.eps')
                plt.close("all")
    