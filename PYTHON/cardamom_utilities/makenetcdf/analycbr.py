2#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:03:58 2018

@author: gquetin

read cbr
analyze cbr parameters
Ultimately save out as a gridded netcdf

"""

#import scipy.io as sio
import numpy as np
import xarray as xr


import os
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 

import cardamom_map_utility as cardmap
import readwritebinary as rwbin
import readcardamommodel as rcardmdl

imp.reload(rwbin)
imp.reload(cardmap)
imp.reload(rcardmdl)


def selectfile(pt,filelist):
    
    filename = []
    for fl in filelist:
        if pt in fl:
            filename.append(fl)
            
    if len(filename) == 1:
        filename = filename[0]
    
    return filename





def cbr_to_netcdf(filelist,cutfrac=.5,INFO=[],filesize=132000,modelid=811):
    '''
    Take a list of a list of files, and the same number of names of and experiment
    save netcdfs of the shared points (can be run with one list)
    Improvements: no need of cheat list, handling multiple runs _1,_2
    '''
    
    
    ptsincommon = list(set([s.split(".")[0].split("_")[-2] for s in filelist]))
    pointlist = np.array([s.split(".")[0].split("_")[-2] for s in filelist])
    
    chainsetcounts = []
    for pt in ptsincommon:
        
        chainsetcounts.append(np.sum(pointlist == pt))
        
    maxchaincount = np.max(chainsetcounts)   
    
    ensemsize = 250
    ensemnum = ensemsize * maxchaincount
    
    
    #
    # Create coordinates to include whole set of points        
    lat_list = []
    lon_list = []
    for colrownum in ptsincommon:
       
        latlontup = cardmap.rowcol_to_latlon([colrownum])
        
        lat, lon = latlontup[0]
        lat_list.append(lat)
        lon_list.append(lon)
        
    latstep = 4
    lonstep = 5
    
    latarray = np.array(lat_list)
    latcoord = np.arange(latarray.min(),latarray.max()+latstep,latstep)
    lonarray = np.array(lon_list)
    loncoord = np.arange(lonarray.min(),lonarray.max()+lonstep,lonstep)
#    ensemnum = 1000
    ensemcord = np.arange(1,ensemnum+1,1)
    
    
    nullmat = np.ones((latcoord.shape[0],loncoord.shape[0],ensemnum))*np.nan
    
    da = xr.DataArray(data=nullmat.copy(),coords = [latcoord,loncoord,ensemcord],
                      dims=['lat','lon','ensemble'])
    

    ds = xr.Dataset()
    
    
    for pn in rcardmdl.get_parnames(modelid,output='shortnames'):
        ds[pn] = da.copy()
        
        
    

    #%

    for cw in ptsincommon:
        
        latlontup = cardmap.rowcol_to_latlon([cw])
        
        lat, lon = latlontup[0]
        
        cbrtemp = selectfile(cw,filelist)

        if type(cbrtemp) == str:
            cbr=[cbrtemp]
        elif type(cbrtemp)==list:
            cbr = sorted(cbrtemp)
        else:
            continue
            print('no cbr here?')
        
        
        
        
        if type(cbr) == list:
            
            pars = []
            for cb in cbr:
                
                if os.path.exists(cb) and os.path.getsize(cb) == filesize:
                    parsraw = rwbin.read_cbr_file(cb,INFO)
                    pars.append(parsraw[int(cutfrac*parsraw.shape[0]):,:])
                    
#                else: #(removed because of pre-filtering)
#                    pars.append([])
            
            if pars == []:
                continue
            
#            #
#            ## Fill in too small files (removed because of pre-filtering)
#            for pr in pars:
#                if not len(pr) == 0:
#                    nshape = pr.shape
#            
#            for jjy, pr in enumerate(pars):
#                if len(pr) == 0:
#                    pars[jjy] = np.ones(nshape)*np.nan
            
            
            pars = np.concatenate(pars,axis=0)
            
            if pars.shape[0] < ensemnum:
                
                newshape = (ensemnum - pars.shape[0],pars.shape[1])
                nullpars = np.ones(newshape)*np.nan
                
                pars = np.concatenate([pars,nullpars],axis=0)
                
        else:
        
            parsraw = rwbin.read_cbr_file(cbr,INFO)
            pars = parsraw[:int(cutfrac*parsraw.shape[0]),:]
    
        for pn in rcardmdl.get_parnames(modelid,output='shortnames'):
            ds[pn].loc[dict(lat=lat,lon=lon)] = np.array(pars[:,rcardmdl.get_parnames(modelid,output='shortnames').index(pn)])
            ds[pn].attrs = {'longname':rcardmdl.get_parnames(modelid=modelid,output='longnames')
                                        [rcardmdl.get_parnames(modelid,output='shortnames').index(pn)]}#,
                              #'units':parnames(modelid=modelid,output='dictunits')[fl]}

    return ds


        
        
def dataset_to_array(ds,varlist):
    datalist = []
    for pn in varlist:
        
        
        datalist.append(ds[pn].values)
        
    dataarray = np.stack(datalist,axis=1)
    
    return dataarray



def netcdf_to_cbr(data_set,nametag,modelid,xy_name = ['lat','lon'],num='0'):
    
    var_list = rcardmdl.get_parnames(modelid,output='shortnames')
    parall_array = np.stack([data_set[var].values for var in var_list]).T
    
    if np.all(~np.isnan(np.nanmean(parall_array,axis=0))):
        
        # Remove nan parameter sets
        param_goodvec = ~np.isnan(parall_array.mean(axis=1))
        
        lat = data_set[xy_name[0]].values
        lon = data_set[xy_name[1]].values
        ptname = cardmap.latlon_to_rowcol(lat,lon)
        
        
        outputfilename = '_'.join([nametag,ptname,num])+'.cbr'
        
        rwbin.write_cbr_file(parall_array[param_goodvec,:], outputfilename)
        success = 1
    else:
        success = 0
        
        
    return success


def netcdfmap_to_cbrs(data_set_xy,nametag,modelid = 811):
    
    lons,lats = np.meshgrid(data_set_xy['lon'].values,data_set_xy['lat'].values)
    
    yesno = []
    for lat,lon in zip(lats.flatten(),lons.flatten()):
        
        data_set = data_set_xy.sel(lat=lat,lon=lon)
        
        yesno.append(netcdf_to_cbr(data_set,nametag,modelid))
        
    return yesno


  
# %%

if __name__ == '__main__':
    
    print('No main file functions currently, check cbrdir2netcdfmap.py')