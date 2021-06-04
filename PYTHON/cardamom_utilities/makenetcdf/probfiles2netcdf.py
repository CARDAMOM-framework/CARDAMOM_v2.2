#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 17:46:21 2019

@author: gquetin


Collect probabilities for a CARDAMOM run into a map

"""

import numpy as np
import xarray as xr

import os
import glob
import sys
import imp

MACHDIR = os.path.expanduser("~") 
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5')

import cardamom_map_utility as cardmap
import file_organization_utility as forguti
import readwritebinary as rwbin

imp.reload(rwbin)
imp.reload(cardmap)
imp.reload(forguti)

def makeblank_geoschem(pts,dictset = None):

    # Default
    if dictset is None:
        ensemnum = 2500
        latstep = 4
        lonstep = 5
        dims = ['lat','lon','ensemble']
    else:
        ensemnum = dictset['ensemnum']
        latstep = dictset['latstep']
        lonstep = dictset['lonstep']
        dims = dictset['dims']

    # Create latitudes to span the points
    lat_list,lon_list = zip(*cardmap.rowcol_to_latlon(pts))
    latarray = np.array(lat_list)
    latcoord = np.arange(latarray.min(),latarray.max()+latstep,latstep)
    lonarray = np.array(lon_list)
    loncoord = np.arange(lonarray.min(),lonarray.max()+lonstep,lonstep)
    
    #
    ## Create a 'probability' or 'ensemble' member vector
    ensemcord = np.arange(1,ensemnum+1,1)


    #
    ## Create a dataarray of nan values to be filled
    nullmat = np.ones((latcoord.shape[0],loncoord.shape[0],ensemnum))*np.nan
    
    da = xr.DataArray(data=nullmat.copy(),coords = [latcoord,loncoord,ensemcord],
                      dims=dims)
    
    return da

def filename_to_unique_colrow(filelist):
    pts = list(set([s.split(".")[0].split("_")[-2] for s in filelist]))
    
    return pts

def filename_to_wildcard(filelist, uniquename = None):
    
    wildcard = list(set(["_".join(s.split("/")[-1].split(".")[0].split("_")[1:-2]) for s in filelist]))
    
    if uniquename is not None and len(wildcard) > 1:
        wildcard_out = [wc for wc in wildcard if uniquename in wc][0]
    elif len(wildcard) > 1:
        print('Too many file names, provide unique name selection')
        wildcard_out = None
    else:
        wildcard_out = wildcard[0]
        
    return wildcard_out

def filename_to_directory(filelist, uniquename = None):
    
    dirname = list(set(["/".join(s.split("/")[:-1]) for s in filelist]))
    
    if uniquename is not None and len(dirname) > 1:
        dirname_out = [dn for dn in dirname if uniquename in dn][0]
    elif len(dirname) > 1:
        print('Too many file names, provide unique name selection')
        dirname_out = None
    else:
        dirname_out = dirname[0]
        
    return dirname_out

def probfiles_to_netcdf(filelist,cutfrac=.5):
    '''
    Take a list of a list of files, and the same number of names of and experiment
    save netcdfs of the shared points (can be run with one list)
    Improvements: no need of cheat list, handling multiple runs _1,_2
    '''
    varnotime = 'probs'
    pfnamestart = 'probfile'
    pts = filename_to_unique_colrow(filelist)
    wildcard = filename_to_wildcard(filelist)
    dirused = filename_to_directory(filelist)
    da =  makeblank_geoschem(pts)
    
    ds = xr.Dataset()
    ds[varnotime] = da.copy()
    
    
    #
    ## Loop through the point list
    for cw in pts:
        
        latlontup = cardmap.rowcol_to_latlon([cw])
        
        lat, lon = latlontup[0]
        
        probfilename = '_'.join([pfnamestart,wildcard,cw])
        probfiles = sorted(glob.glob(dirused + '/' + probfilename + '*.bin'))
        
        if len(probfiles) == 0:
            print('no probfile here?')
            continue
            
        probs = []
        for pb in probfiles:
            
            if os.path.exists(pb):
                probraw = rwbin.readbinarymat(pb,[1])
                probs.append(probraw[int(cutfrac*probraw.shape[0]):])
        
        probs = np.concatenate(probs,axis=0)
        
        ds[varnotime].loc[dict(lat=lat,lon=lon,ensemble=slice(1,np.array(probs).shape[0]))] = np.array(probs).squeeze()
        ds[varnotime].attrs = {'longname':'probfile'}
    
        

    return ds


# %%

if __name__ == '__main__':
    
    # Set up folders based on the home folder
    datadir, dataoutdir = forguti.folderdeff(MACHDIR)
    
    exp = 'p20_long' #'p12_uncmatch_full' #'p16_long'
    driver = 'cruncep' #'cruncep_osse5'
    diroutput = '{}DALEC/{}/{}/output/'.format(datadir,exp,driver)
#    diroutput = '/Users/gquetin/Google Drive/CARDAMOM_Share/SharedNotes_Experiments/ball_berry_model_tests/data/probfiles/'
    
    # Test the creation of a probfile
    filelist = glob.glob(diroutput + 'probfile*.bin')
    dstest = probfiles_to_netcdf(filelist,cutfrac=.5)
    
    dstest.to_netcdf(diroutput + 'SUMMARYSTAT1_probfile_cru004GCR006_1920_2015_nbe2002.nc')