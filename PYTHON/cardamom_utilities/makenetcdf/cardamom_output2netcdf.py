#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:37:39 2019

@author: gquetin

Load flux and pool files of CARDAMOM output to netcdf files

Options:
    Filtering:
        File size
        Zero variance
    Variables (split pools and fluxes?)
    How to deal with partial files
    Make it submittable to Sherlock in parallel [information from commandline]

"""

import scipy.io as sio
import numpy as np
import xarray as xr




import os
import glob
import sys
import imp

MACHDIR = os.path.expanduser("~") 
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5')

import file_organization_utility as forguti
import cardamom_map_utility as cardmap
import readwritebinary as rwbin
import makecbf


imp.reload(forguti)
imp.reload(cardmap)
imp.reload(rwbin)
#import dalecmetfileload as dmet

#import dalecmetfileload as dmet

import outputnames
    

def cbrmat_to_netcdf(filename,lat,lon,time,baddata=-9999):
                
    """
    Load .mat or .bin file into a netcdf for easier reading/storing
    """
    
    # DEFINITIONS
    pools = outputnames.poolnames(output='shortnames')

    fluxes = outputnames.fluxnames(output='shortnames')

    # Load .mat file output from DALEC
    if type(filename) == str and filename.split(".")[-1] == 'mat':
        dalecarray = sio.loadmat(filename,struct_as_record=False, squeeze_me=True)
        output_struct = dalecarray['CBR']
        POOLS = output_struct.POOLS
        FLUXES = output_struct.FLUXES
        LAI = output_struct.LAI
        GPP = output_struct.GPP
        NBE = output_struct.NBE
        NEE = output_struct.NEE
        
    elif type(filename) == dict:
        POOLS = filename['POOLS']
        FLUXES = filename['FLUXES']
        
        if 'PROBS' in list(filename.keys()):
            PROBS = filename['PROBS']
        LAI = filename['LAI']
        GPP = filename['GPP']
        NBE = filename['NBE']
        NEE = filename['NEE']
        
    else:
        print('Not recognized format')
        return
        
    #
    ## Clean up bad data
    POOLS[POOLS == baddata] = np.nan
    FLUXES[FLUXES == baddata] = np.nan
    if 'PROBS' in list(filename.keys()):
        PROBS[PROBS == baddata] = np.nan
    LAI[LAI == baddata] = np.nan
    GPP[GPP == baddata] = np.nan
    NBE[NBE == baddata] = np.nan
    NEE[NEE == baddata] = np.nan


    # Vector of parameters sets, just labeled 1 to 1000
    diffpars = np.arange(1,POOLS.shape[0]+1,1) #number of models run with differen parameters
    


    #
    # Create dataset of outputs
    dspools = xr.Dataset()
    
    for idx, pl in enumerate(pools):
        dspools[pl] = xr.DataArray(POOLS[:,:,idx],coords=[diffpars,time],
               dims = ['parameterset','time'])
        dspools[pl].attrs = {'longname':outputnames.poolnames(output='dict')[pl],
                             'units':outputnames.poolnames(output='dictunits')[pl]}
    
    
    
    dsfluxes = xr.Dataset()
    for idx, fl in enumerate(fluxes):
        dsfluxes[fl] = xr.DataArray(FLUXES[:,:,idx],
                coords={'parameterset':diffpars,'time':time,
                        'lat':lat,'lon':lon},
                dims = ['parameterset','time'])
        dsfluxes[fl].attrs = {'longname':outputnames.fluxnames(output='dict')[fl],
                              'units':outputnames.fluxnames(output='dictunits')[fl]}
    if 'PROBS' in list(filename.keys()):
            
        daprobs = xr.DataArray(PROBS.squeeze(),
                               coords={'parameterset':diffpars,
                                       'lat':lat,'lon':lon},
                                       dims=['parameterset'],name='probs')
    
    
    #    dafluxes = xr.DataArray(output_struct.FLUXES,coords=[diffpars,time,fluxes],
    #                            dims=['parameterset','time','flux'],name='fluxes')
    dalai = xr.DataArray(LAI,
                         coords={'parameterset':diffpars,'time':time,
                        'lat':lat,'lon':lon},
                         dims=['parameterset','time'],name='lai')
    dalai.attrs['longname'] = outputnames.derivedfluxnames(output='dict')['lai']
    dalai.attrs['units'] = outputnames.derivedfluxnames(output='units')
    
    dagpp = xr.DataArray(GPP,
                         coords={'parameterset':diffpars,'time':time,
                        'lat':lat,'lon':lon},
                         dims=['parameterset','time'],name='gpp')
    dagpp.attrs['longname'] = outputnames.derivedfluxnames(output='dict')['gpp']
    dagpp.attrs['units'] = outputnames.derivedfluxnames(output='units')
    
    danbe = xr.DataArray(NBE,
                         coords={'parameterset':diffpars,'time':time,
                        'lat':lat,'lon':lon},
                         dims=['parameterset','time'],name='nbe')
    danbe.attrs['longname'] = outputnames.derivedfluxnames(output='dict')['nbe']
    danbe.attrs['units'] = outputnames.derivedfluxnames(output='units')
    
    danee = xr.DataArray(NEE,
                         coords={'parameterset':diffpars,'time':time,
                        'lat':lat,'lon':lon},
                         dims=['parameterset','time'],name='nee')
    danee.attrs['longname'] = outputnames.derivedfluxnames(output='dict')['nee']
    danee.attrs['units'] = outputnames.derivedfluxnames(output='units')
    
    if 'PROBS' in list(filename.keys()):   
        dsdalecoutput = xr.merge([dspools, dsfluxes,dalai,dagpp,danbe,danee,daprobs])
    else:
        dsdalecoutput = xr.merge([dspools, dsfluxes,dalai,dagpp,danbe,danee])
    
    return dsdalecoutput



def checkbadparm(dataset):

    dscheck = dataset.mean(dim='time',skipna=False)

    badparm = []
    for vc in np.array(dscheck.data_vars):
        
        # Loop through and identify any parameters with an infinity
        
        badparm.append(dscheck['parameterset'].values[np.abs(dscheck[vc].values) > 1e6])
        badparm.append(dscheck['parameterset'].values[np.isnan(dscheck[vc].values)])
    
    badparmlist = np.array(set(np.concatenate(badparm))).tolist()
    goodparmlist = [s for s in dscheck['parameterset'].values if s not in badparmlist]
    
    datasetgood = dataset.sel(parameterset=goodparmlist)
    
    return datasetgood

#fluxes = ['Flux'+ str(s) for s in np.arange(1,output_struct.FLUXES.shape[2]+1,1)]

#POOLS = {'Woody Biomass':3,'Plant Available Water':6}

def file_to_colrow(filename,call = None):
    '''
    Converts standard CARDAMOM output name to the 
    column+row number, or the wildcard name (everything else)
    
    <typefile>_<wildcard>_<colrow>_<chain>.<suffix>
    
    typefile = poisition 0
    suffix = after '.'
    colrow = position -2 if suffix = .cbr, .bin or position -1 if .cbf
    chain = position -1 if suffix = .cbr, .bin
    '''
    
    name1 = filename.split("/")[-1].split(".")
    if name1[-1] == 'bin':
        colrownum = name1[0].split("_")[-2]
        wildcard = "_".join(name1[0].split("_")[1:-2])
        chain = name1[0].split("_")[-1]
    
    elif name1[-1] == 'cbr':
        colrownum = name1[0].split("_")[-2]
        wildcard = "_".join(name1[0].split("_")[0:-2])
        chain = name1[0].split("_")[-1]
            
    elif name1[-1] == 'cbf':
        colrownum = name1[0].split("_")[-1]
        wildcard = "_".join(name1[0].split("_")[0:-1])
    
    if call == 'colrow':
        return colrownum
    elif call == 'wildcard':
        return wildcard
    elif call == 'chain':
        return chain
    else:
        return colrownum, wildcard





def checklistoverlap(checklist):
    count = 0
    while count < len(checklist)-1:
        if count == 0:
            setcheck = checklist[0]
        setcheck = list(set(setcheck) & set(checklist[count+1]))
        count = count + 1
        
    return setcheck



def filedict_to_netcdf(filedict,sc,cutfrac=0.5):
    latlontup = cardmap.rowcol_to_latlon([sc])
    lat, lon = latlontup[0]
    
    typelist = [tyl for tyl in list(filedict.keys()) if tyl not in ['cbr','cbf']]
    fbins = []
    for ky in typelist:
        
        fbins += [fl for fl in filedict[ky] if file_to_colrow(fl,call='colrow') == sc]
    
    #
    ## Check for run numbers
    runnumlist = list(set([fb.split("/")[-1].split("_")[-1].split(".")[0] for fb in fbins]))
    
    fcbrsall = [fl for fl in filedict['cbr'] if file_to_colrow(fl,call='colrow') == sc]
    fcbrs = [fb for fb in fcbrsall if fb.split("/")[-1].split("_")[-1].split(".")[0] in runnumlist]
    
    if len(runnumlist) != len(fcbrs):
        print('Error: point ('+sc+') missing a .cbr file')
        return
    
    fcbf = [fl for fl in filedict['cbf'] if file_to_colrow(fl,call='colrow') == sc]
        
    if len(fcbf) != 1:
        print('Error: point ('+sc+') too many .cbf files')
        return
    
    CBF = rwbin.read_cbf_file(fcbf[0])
    

    CBR = rwbin.CARDAMOM_READ_OUTPUT_FILEAWARE(fcbf[0],fcbrs,fbins,cutfrac = cutfrac)
        
    
    
    metlabels = makecbf.cbflabels()
    dayarray = CBF['MET'][:,metlabels.index('Time [Days since Jan 01, 2001]')]
    time = makecbf.datetime_from_cbf(dayarray,refdate='2001-01-01') 
    
    dsdalec = cbrmat_to_netcdf(CBR,lat,lon,time)
    
    return dsdalec


def time_from_cbf(CBF,tbase = '2001-01-01',tbasename = 'Time [Days since Jan 01, 2001]'):
    metlabels = makecbf.cbflabels()
    dayarray = CBF['MET'][:,metlabels.index(tbasename)]
    time = makecbf.datetime_from_cbf(dayarray,refdate=tbase) 
    
    return time


def allfile_ptlist(dirdict, wcdict, typelist = ['probfile','edcdfile','fluxfile','poolfile']):
    '''
    Return set of points that are in all .cbr, .cbf, and .bin files provided in dirdict and wcdict
    '''
    
    filedict = {}
    for typefile in typelist:
        
        try:
            wildcard = wcdict['output']
            filedict[typefile] = sorted([fl for fl in glob.glob(dirdict['output']+typefile+'_'+wildcard + '*.bin')])
        except:
            print('No '+ typefile + ' available')
            
        
        
        
    filedict['cbf'] = sorted(glob.glob(dirdict['cbf']+wcdict['cbf']+'*.cbf'))     
    filedict['cbr'] = sorted(glob.glob(dirdict['cbr']+wcdict['cbr']+'*.cbr'))
    
    
    colrowlist = []
    for ky in filedict:  
        tmplist=[]
        
        for fl in filedict[ky]:
            if type(fl)== str:
                tmplist.append(file_to_colrow(fl,call='colrow'))
            elif type(fl)==list:
                tmplist.append(file_to_colrow(fl[0],call='colrow'))

        colrowlist.append(tmplist)
    
    
    setcheck = checklistoverlap(colrowlist)
    
    return setcheck





def open_fluxes_ptdir(sc, dirdict, wcdict,settings,typelist):

    '''
    Receives:
        sc - the grid point
        dirdict - directories of where the various files are and go
        wcdict - wildcards to identify files where there are multiple in a directory
        settings - which chains at that point to process and what fraction of each chain
        typelist - all the files to process
    
    '''
    
    if 'cutfrac' in settings and settings['cutfrac'] is not None:
        cutfrac = settings['cutfrac']
    else:
        cutfrac = 1.0

    

    fcbf = "_".join([dirdict['cbf']+wcdict['cbf'],sc])+'.cbf'
    
    fcbrlist = []
    livelist = glob.glob(dirdict['cbr']+wcdict['cbr']+'_'+sc+'*.cbr')
    cbrchainset = []
    for fn in livelist:
        
        if 'chains' in list(settings.keys()) and settings['chains'] is not None:
            
            pickchains = settings['chains']
            
            chain = file_to_colrow(fn,call = 'chain')
            if int(chain) in pickchains:
                fcbrlist.append(fn)
                cbrchainset.append(chain)
        
        else:
            chain = file_to_colrow(fn,call = 'chain')
            fcbrlist.append(fn)
            cbrchainset.append(chain)
            
    fcbrlist = sorted(fcbrlist)
            
        
    #
    ## Check that .cbr files exist
    
    fcbrlist_checked = [fcbr for fcbr in fcbrlist if fcbr in livelist]
    cbrchainset = set([file_to_colrow(ch,call = 'chain') for ch in fcbrlist_checked])
    
    fbinlist = []
    binchainset = []
    for typefile in typelist:
        livelist = glob.glob(dirdict['output']+typefile+'_' + wcdict['output']+'_'+sc+'*.bin')
        
        for fn in livelist:
            
            if 'chains' in list(settings.keys()) and settings['chains'] is not None:
            
                pickchains = settings['chains']
                
                chain = file_to_colrow(fn,call = 'chain')
                if int(chain) in pickchains:
                    fbinlist.append(fn)
                    binchainset.append(chain)
        
            else:
                chain = file_to_colrow(fn,call = 'chain')
                fbinlist.append(fn)
                binchainset.append(chain)
                
    fbinlist = sorted(fbinlist)
            
    
    #
    ## Check that .bin files exist
    allchains = list(set(binchainset) & set(cbrchainset)) 
    fcbrlist_checked = [fcbr for fcbr in fcbrlist if file_to_colrow(fcbr,call = 'chain') in allchains]
    fbinlist_checked = [fbin for fbin in fbinlist if file_to_colrow(fbin,call = 'chain') in allchains]

    
    CBF = rwbin.read_cbf_file(fcbf)
    
    CBR, Nlist = rwbin.CARDAMOM_READ_OUTPUT_FILEAWARE(fcbf,fcbrlist_checked,fbinlist_checked,cutfrac = cutfrac,chaininfoout='list')
    
    return CBF, CBR, Nlist


def output2savenetcdf(sc, dirdict, wcdict,settings,typelist,varpick = None, datasave=True):
        
    if not os.path.exists(dirdict['output']+"_".join([wcdict['output'],sc])+'.nc'):

        try:
            
            CBF, CBR, Nlist = open_fluxes_ptdir(sc, dirdict, wcdict,settings,typelist)
            
            lat, lon = cardmap.rowcol_to_latlon([sc])[0]
            time = time_from_cbf(CBF)
            
            dsdalec = cbrmat_to_netcdf(CBR,lat,lon,time)
            dsdalec.attrs['chainnum'] = Nlist
            
            #
            ## Select for fewer variables
            if varpick is None:
                dsdalec_tosave = dsdalec
                
            else:
                dsdalec_tosave = dsdalec[varpick]
            
            
            # Save out to netCDF
            if datasave is True:
                dsdalec_tosave.to_netcdf(path=dirdict['output']+"_".join([wcdict['output'],sc])+'.nc')
                
            return dsdalec_tosave
            
        except:
            print('badpoint '+sc)
            return None
        
    else:
        print('File already exists for ' + sc)
        return None

# %%

if __name__ == '__main__':

    
    datasave = True
    
    # Set up folders based on the home folder
    datadir, dataoutdir = forguti.folderdeff(MACHDIR)
    

    mdl = 'cruncep' #'cesmlensv3bc001'#'ctrllensbc001'#'ctrllensbc005' #'ctrllensbc001' #'cruncep' #'HadGEM2-ES' #'erainterim'#[MODELPICKS[0]]
    mdl_other = 'cruncep' #'cruncep_forward' #'cruncep_forward'
    experm = 'p20_long' #'p16_long' #'p12_uncmatch_full' #'p16_long' #'p12_uncmatch_lens001' #'p3_uncmatch_limit_20' #'multipt_test'     
    experm_other = 'p20_long'

    # Default setting
    settings = {'cutfrac': 0.0, #0.5,
                'chains': [0]}#[1,2,3,4,5,6,7,8,9,10]}
   
    varpick = None
    
    dirdict = {}
    dirdict['output'] = '{}DALEC/{}/{}/{}/'.format(datadir,experm,mdl,'output')
    dirdict['cbf'] = '{}DALEC/{}/{}/{}/'.format(datadir,experm,mdl,'cbf')
    dirdict['cbr'] = '{}DALEC/{}/{}/{}/'.format(datadir,experm_other,mdl_other,'cbr')
    
    wcdict = {}
    wcdict['output'] = file_to_colrow(glob.glob(dirdict['output']+'*.bin')[0],call='wildcard')
    wcdict['cbf'] = file_to_colrow(glob.glob(dirdict['cbf']+'*.cbf')[0],call='wildcard')
    wcdict['cbr'] = file_to_colrow(glob.glob(dirdict['cbr']+'*.cbr')[0],call='wildcard')
    wcdict['probfile'] = wcdict['output']
    

#    typelist = ['probfile','edcdfile','fluxfile','poolfile']
    
    setcheck = allfile_ptlist(dirdict, wcdict)
    
    #%% 
    #
    ## Need a function that loops through a list of point + chain and reads out points
    ## for setcheck and chains to be used
    ## Could put a list of chains on the command line
    
    #
    ## Make wcdict and dirdict as part of the command line input?
    
    #
    ## If looping whole file use this to find good points

    input_list = sys.argv
    print(input_list)
    
    if len(input_list) > 1: # 'pt'
        setcheck = [str(input_list[1])]
        print('Commandline inputs: ' + input_list[1])

        
    
    if len(input_list) > 2: # 'chainlist'
        chainstring = input_list[2]
        chains = [int(jj) for jj in chainstring.split("_")]
        settings['chains'] = chains
        
    # Set the cutfraction for the analysis
    if len(input_list) > 3: # 'cutfrac'
        settings['cutfrac'] = float(input_list[3])
        
    if len(input_list) > 4: # 'varcall' list of variables to save out in the netcdf
        varcall = str(input_list[4])
        
        fluxpick = outputnames.derivedfluxnames(output='shortnames') + ['et','gpp_to_autoresp','hetresp_litter','hetresp_som','fire_em_total']
        poolpick = outputnames.poolnames(output='shortnames')
        
        var_dict = {'poolflux': fluxpick + poolpick,
                    'flux':fluxpick,
                    'pool':poolpick}
        
        try:
            varpick = var_dict[varcall]
            
        except:
            print('Not an option: ' + varcall)
            varpick = None
        
    if len(input_list) > 5: # 'cbfdict'
        dirdict['cbf'] = str(input_list[5])
        wcdict['cbf'] = file_to_colrow(glob.glob(dirdict['cbf']+'*.cbf')[0],call='wildcard')
    
    if len(input_list) > 6: # 'cbrdict'
        dirdict['cbr'] = str(input_list[6])
        wcdict['cbr'] = file_to_colrow(glob.glob(dirdict['cbr']+'*.cbr')[0],call='wildcard')
        
    if len(input_list) > 7: # 'outdict'
        dirdict['output'] = str(input_list[7])
        wcdict['output'] = file_to_colrow(glob.glob(dirdict['output']+'*.bin')[0],call='wildcard')
        wcdict['probfile'] = wcdict['output']
        
    
    # Take out the probifle if calculating for forward run only
    if wcdict['cbf'] != wcdict['cbr']:
        typelist = ['edcdfile','fluxfile','poolfile']
    else:
        typelist = ['probfile','edcdfile','fluxfile','poolfile']
    

    
    
    
    

    

    
    #%% Run the processing script
    # Input can be any lenght string of points
    # Make settings iterable per point
    
    #
    # Calcs, create point netCDFs and a netCDF of calculations
    for sc in setcheck:
        print(sc)
        dsdalec = output2savenetcdf(sc, dirdict, wcdict,settings,typelist,varpick=varpick, datasave=datasave)
        
        
        
                    
                        
                    
                    
            



        
