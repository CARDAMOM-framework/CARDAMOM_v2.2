#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:38:52 2020

@author: gquetin

Make a summary plot to check .cbf file input
"""

import readwritebinary as rwbin
import makecbf
import cardamom_output2netcdf as carop2nc

import numpy as np
import os
import glob
import sys

import matplotlib.pyplot as plt

MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5')






def loadcbf(inputfilename = 'test'):
    #
    ## Test file
    if inputfilename == 'test':
        dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p20_long/cruncep/cbf/'
        fnbinary = 'cru004GCR006_1920_2015_nbe2002_2042.cbf'
        inputfilename = dirbinary+fnbinary
    
    
    CBF = rwbin.read_cbf_file(inputfilename)
    
    return CBF





#%%

def checkplotmet(inputfilename,plotlocation = None):
    CBF = loadcbf(inputfilename = inputfilename)
    time = carop2nc.time_from_cbf(CBF)
    MET = CBF['MET']
    metlabels = makecbf.cbflabels()
    
    fig, axarr = plt.subplots(MET.shape[1],1,figsize=(10,10),sharex=True)
    
    for jjm, met in enumerate(MET.T):
        
        ax = axarr.flatten()[jjm]
        
        ax.plot(time,met)
        if jjm == 0:
            ax.set_title('{}\n{}'.format(inputfilename.split('/')[-1],metlabels[jjm]))
        else:
            ax.set_title(metlabels[jjm])
        
        if jjm == MET.shape[1]-1:
            ax.set_xlabel('time')
        
    fig.tight_layout()
    plt.show()
    
    return fig


#%%

def checkplotobs(CBF,plotlocation = None):
    time = carop2nc.time_from_cbf(CBF)
    OBSdict = CBF['OBS']
    fig, axarr = plt.subplots(len(OBSdict),1,figsize=(6,10),sharex=True)
    
    for jjm, kyob in enumerate(OBSdict):
        
        ax = axarr.flatten()[jjm]
        
        data = np.array(OBSdict[kyob])
        
        if kyob == 'GPP' and CBF['OBSUNC']['GPP']['gppabs'] == 0.0:
            ax.set_title('SIF/Proportional')
        else:
            ax.set_title(kyob)
        
        if data.shape[0] != time.shape[0]:
            print('{} not used'.format(kyob))
            
            continue
        
        
        
        data[data==-9999.] = np.nan
        
        ax.plot(time,data)
        
        
        
        if jjm == len(OBSdict)-1:
            ax.set_xlabel('time')
        
    fig.tight_layout()
    plt.show()
    return fig


def checkall(inputfilename='test',plotlocation=None):

    CBF = loadcbf(inputfilename)
    
    # Plot met
    checkplotmet(CBF,plotlocation = None)
    
    # Plot obs
    checkplotobs(CBF,plotlocation = None)
    

# %%

if __name__ == '__main__':
    # args?
    EXP = 'p20_long'
    dirtop = '{}/Google Drive/DATA/DALEC/{}/'.format(MACHDIR,EXP)
    ptcount_test = 20 #944
    
    direxp_list = glob.glob(dirtop+'*')
    mettype = 'ctlensv3bc001'#'upco2ctlv3bc' #'noco2lnsv3bc' #'cesmlensv3bc'
    
    direxp_list = [dx for dx in direxp_list if mettype in dx]
    obcheck = 'LAI' #'ABGB'
    for drex in direxp_list[0:1]:
        
        cbf_filelist = glob.glob('{}/cbf/*.cbf'.format(drex))
        
        
        cbf_filelist_totest = np.random.permutation(cbf_filelist)[:ptcount_test]
    
        for inputfilename in cbf_filelist_totest:
#            print(inputfilename.split('/')[-1])
            
            CBF = loadcbf(inputfilename)
            
            A = CBF['OBS'][obcheck].copy()
            if A != []:
                A[A == -9999] = np.nan
                print('{}: {}'.format(obcheck,np.nanmean(A)))
            checkplotmet(inputfilename)
#            checkplotobs(CBF)
#            checkall(inputfilename=inputfilename)