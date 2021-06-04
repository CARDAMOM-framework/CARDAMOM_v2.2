#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:37:34 2018

@author: gregoryquetin


Reading/Writing .cbf binary file
"""

import numpy as np
import scipy.io as sio
import os
import glob
import sys
import warnings

MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5')


import readcardamommodel



 



#% Functions

def readbinarymat(filename,indims):
    
    alldim = 1
    for dm in indims:
        alldim = alldim*dm
        
    with open(filename, 'rb') as fid:
        BD = np.fromfile(fid, np.float64)

    N = int(BD.shape[0]/alldim)

    dims = [N] + indims
    
    try:
        datamat = BD.reshape(dims)
    except:
        print('Not a complete file')
        datamat = None
    
    return datamat



def define_cbf_obs_fields(CBF,CBFOBS):
    '''
    #Time-resolved observations
    '''

    obsnames=['GPP','LAI','NBE','ABGB','ET','EWT','BAND1','BAND2','BAND3','BAND4','SOM']
    
    if 'OBS' in CBF.keys():
        print('Waring: OBS already in CBF')
    else:
        CBF['OBS'] = {}

    
    
    for n, obn in enumerate(obsnames):  #CBF.noobs;
           if np.shape(CBFOBS)[1]>=n+1 and ~np.all(CBFOBS[:,n]== -9999):
               CBF['OBS'][obn]=CBFOBS[:,n]
           else:
               CBF['OBS'][obn]=[];

        
    CBF['OBSinfo'] = {}
    CBF['OBSinfo']['LAI']=('LAI data requirements: should be +ve definite (>0);'+
                           'assumed uncertainty for LAI timeseries is factor of 2;'+
                           'future versions will include user-defined uncertainty options')
    CBF['OBSinfo']['uncertainty_factors']=('Uncertainty structures for positive-definite'+
                                           'quantities (GPP, LAI, ET), are prescribed as '+
                                           'uncertainty factors (by default); uncertainty '+
                                           'factors should be > 1. \n For example: a 1-sigma'+
                                           'range for 100 uncertainty factor 2 = 100/2 - 100*2 = 50 - 200 ')

    return CBF

#Uncertainties for time-resolved observations
def read_obs_uncertainty_fields(CBF,SD,OPRU):
    '''
    %From CARDAMOM_READ_BINARY_DATA.c 
    % DATA->nee_annual_unc=statdat[13];
    % DATA->et_annual_unc=statdat[14];
    % DATA->nee_obs_unc=statdat[15];if (statdat[15]<0){DATA->nee_obs_unc=0.5;}
    % DATA->et_obs_unc=statdat[16];if (statdat[16]<0){DATA->et_obs_unc=2;}
    % DATA->ewt_annual_unc=statdat[17];
    % DATA->ewt_obs_unc=statdat[18];if (statdat[18]<0){DATA->ewt_obs_unc=50;}
    % DATA->gpp_annual_unc=statdat[19];
    % DATA->gpp_obs_unc=statdat[20];if (statdat[18]<0){DATA->gpp_obs_unc=2;}
    '''
    
    if 'OBSUNC' in CBF.keys():
        print('Waring: OBSUNC already in CBF')
    else:
        CBF['OBSUNC'] = {}
    
    #NBE
    CBF['OBSUNC']['NBE'] = {}
    CBF['OBSUNC']['NBE']['annual_unc']=SD[13];
    CBF['OBSUNC']['NBE']['seasonal_unc']=SD[15];
    CBF['OBSUNC']['NBE']['info']='Single point (default = 0.5, must be >0) and annual (annual_unc) NBE uncertainty [gC/m2/d]';
    
    #ET 
    CBF['OBSUNC']['ET'] = {}
    CBF['OBSUNC']['ET']['annual_unc']=SD[14];
    CBF['OBSUNC']['ET']['unc']=SD[16];
    CBF['OBSUNC']['ET']['obs_unc_threshold']=SD[21] #Default = 0.1
    CBF['OBSUNC']['ET']['info']=('Single point uncertainty factor (default = */ +2, must be >1) and annual (annual_unc) ET uncertainty [factor]. \nDefault obs threshold = 0.1 mm/day: this ensures log-tranformed model ET values are insensitive to ET<0.1');
    
    #EWT
    CBF['OBSUNC']['EWT'] = {}
    CBF['OBSUNC']['EWT']['annual_unc']=SD[17];
    CBF['OBSUNC']['EWT']['unc']=SD[18];
    CBF['OBSUNC']['EWT']['info']=('Single point (default = 50) and annual (N/A yet) EWT uncertainty [mm]')
    
    #GPP
    CBF['OBSUNC']['GPP'] = {}
    CBF['OBSUNC']['GPP']['annual_unc']=SD[19];
    CBF['OBSUNC']['GPP']['unc']=SD[20];
    CBF['OBSUNC']['GPP']['info']=('Single point uncertaint factor (default = */ +2, must be >1) and annual (annual_unc) GPP uncertainty [gC/m2/d]\n default obs unc threshold is 0.1gC/m2/d: this ensures log-tranformed model GPP values are insensitive to GPP<0.1');
    #gpp abs: option on treating GPP constraint as absolute or relative
    CBF['OBSUNC']['GPP']['gppabs']=SD[7];
    CBF['OBSUNC']['GPP']['obs_unc_threshold']=SD[22]   #Default = 0.1
    CBF['OBSUNC']['GPP']['gppabs_info']=('Set to "1" for GPP data, set to "0" for SIF data"')
    
    
    #Time-resolved SOM uncertainty
    CBF['OBSUNC']['SOM'] = {}
    CBF['OBSUNC']['SOM']['unc']=OPRU[7];
    CBF['OBSUNC']['SOM']['info']=('Single point uncertainty on time-resolved SOM, CBF.OBS.SOM')
    
    #Time-resolved biomass uncertainty
    CBF['OBSUNC']['ABGB'] = {}
    CBF['OBSUNC']['ABGB']['unc']=OPRU[1];
    CBF['OBSUNC']['ABGB']['info']=('Uncertainty on time-resolved biomass CBF.OBS.ABGB')
    
    return CBF


 
def read_other_obs_constraints(CBF,OPR,OPRU):
    '''
    other (time-invariant) constraints 
    '''
    
    if 'OTHER_OBS' in CBF.keys():
        print('Waring: OTHER_OBS already in CBF')
    else:
        CBF['OTHER_OBS'] = {}
    #Mean biomass
    CBF['OTHER_OBS']['MBiomass'] = {}
    CBF['OTHER_OBS']['MBiomass']['mean']=OPR[0];
    CBF['OTHER_OBS']['MBiomass']['unc']=OPRU[0];
    CBF['OTHER_OBS']['MBiomass']['info']=('Mean total (above-& below- ground) biomass at time t=0\n'+
                                           'Uncertainty distribution is log-normal distribution')
    
    #Mean fire emissions
    CBF['OTHER_OBS']['MFire'] = {}
    CBF['OTHER_OBS']['MFire']['mean']=OPR[2];
    CBF['OTHER_OBS']['MFire']['unc']=OPRU[2];
    CBF['OTHER_OBS']['MFire']['info']=('Mean total (above-& below- ground) fire emissions '+
                              'for whole simulation \n Use +ve value for uncertainty '+
                              'factor  (log-gaussian distribution) & -ve value for '+
                              'absolute uncertainty (naussian distribution)')
    #Mean LAI
    CBF['OTHER_OBS']['MLAI'] = {}
    CBF['OTHER_OBS']['MLAI']['mean']=OPR[4];
    CBF['OTHER_OBS']['MLAI']['unc']=OPRU[4];
    CBF['OTHER_OBS']['MLAI']['info']=('Mean total LAI for whole simulation \n '+
                                      'Use +ve value for uncertainty factor  '+
                                      '(log-gaussian distribution) & -ve value '+
                                      'for absolute uncertainty (gaussian distribution)\n'+
                                      'If any CBF.OBS.LAI values are presribed, then mean LAI'+
                                      'will be used to constrain these timesteps only.')
    #Mean GPP
    CBF['OTHER_OBS']['MGPP'] = {}
    CBF['OTHER_OBS']['MGPP']['mean']=OPR[5];
    CBF['OTHER_OBS']['MGPP']['unc']=OPRU[5];
    CBF['OTHER_OBS']['MGPP']['info']=('Mean total (above-& below- ground) GPP '+
                                      'for whole simulation \n Use +ve value for '+
                                      'uncertainty factor  (log-gaussian distribution) '+
                                      '& -ve value for absolute uncertainty (naussian distribution)')

    return CBF










def read_cbf_file(inputfilename=None):
    """
    Adapted from A. Bloom matlab script from 02/11/2018
    %.cbf = (c)ARDAMOM (b)INARY (f)ILE
    %.cbr = (c)ARDAMOM (b)INARY (r)ESULTS
    %See also:
        %"CARDAMOM_WRITE_BINARY_FILEFORMAT.m"
        %s
        %Last adapted froo A.A. Bloom 2019/07/24
    """
    
    if inputfilename is None:
        BD = np.ones(23580)*-9999.0
        BD[2] = 1164# nodays
        BD[3] = 9# no of met
        BD[4] = 11# no of obs
    
    else:
        with open(inputfilename, 'rb') as fid:
            BD = np.fromfile(fid, np.float64)
        
    # https://www.mathworks.com/help/matlab/ref/fwrite.html
    # https://www.mathworks.com/help/matlab/numeric-types.html
        
    k=0;
    # Static data (100 places)
    SD=BD[k:k+100]
    k=k+100
    # Priors (50 places)
    PR=BD[k:k+50];
    k=k+50;
    # Priorunc (50 places)
    PRU=BD[k:k+50]
    k=k+50    
    
    # O. Priors (50 places)
    OPR=BD[k:k+50]
    k=k+50
    # O. Priorunc (50 places)
    OPRU=BD[k:k+50]
    k=k+50
    
    CBF = {}
    CBF['PARPRIORS'] = np.expand_dims(PR,axis=1)
    CBF['PARPRIORUNC'] = np.expand_dims(PRU,axis=1)
    CBF=read_other_obs_constraints(CBF,OPR,OPRU)
    
    CBF['ID'] = SD[0] # ID (not used)
    CBF['LAT'] = SD[1] # Latitude
    CBF['nodays'] = int(SD[2]) # Number of days
    CBF['nomet'] = int(SD[3])
    CBF['noobs'] =int(SD[4])
    CBF['EDC'] = SD[5]
    CBF['EDCDIAG'] = SD[6]
#    CBF = {'PARPRIORS':np.expand_dims(PR,axis=1), 
#           'PARPRIORUNC':np.expand_dims(PRU,axis=1), 
#           'OTHERPRIORS':np.expand_dims(OPR,axis=1), #
#           'OTHERPRIORSUNC':np.expand_dims(OPRU,axis=1),
#           'ID':SD[0], # ID (not used)
#           'LAT':SD[1], # Latitude
#           'nodays':int(SD[2]), # Number of days
#           'nomet':int(SD[3]), 
#           'noobs':int(SD[4]),
#           'EDC':SD[5],
#           'EDCDIAG':SD[6],
#           'gppabs':SD[7],
#           'rc_random_search':SD[10]==1,
#           'nbe_annual_unc':SD[13],
#           'etiav':SD[14],
#           'nbe_seasonal_unc':SD[15]}
    
    #MCMC start searching EDCs from anywhere (1) or from prescribed starting
    #point(0). this is outdated - consider deleting
    CBF['rc_random_search'] = SD[10]==1
    
    #NEE IAV options
    CBF=read_obs_uncertainty_fields(CBF,SD,OPRU)
    
    
    TEMPDATA=BD[k:k+(CBF['nomet']+CBF['noobs'])*CBF['nodays']].reshape(CBF['nodays'],
                                                       (CBF['nomet']+CBF['noobs']))
    #All met data
    CBF['MET'] = TEMPDATA[0:CBF['nodays'],0:CBF['nomet']] # Add in new meteorology here
#    CBF['OBS'] = TEMPDATA[0:CBF['nodays'],CBF['nomet']:]
    CBFOBS = TEMPDATA[0:CBF['nodays'],CBF['nomet']:]
    CBF=define_cbf_obs_fields(CBF,CBFOBS)
    
    #Removing redundant fields
#    CBF=rmfield(CBF,'noobs');
#    # CBF=rmfield(CBF,'nomet');
#    # CBF=rmfield(CBF,'nodays');
    
    
    # Read prescribed mean meteorology
    
    if len(BD) - (k+(CBF['nomet']+CBF['noobs'])*CBF['nodays']) == CBF['nomet'] + CBF['noobs']:
        
        kmmet= k+(CBF['nomet']+CBF['noobs'])*CBF['nodays']
        CBF['mmet'] = BD[kmmet:kmmet+CBF['nomet']]
    
    
    #Retaining "OTHERPRIORS" for now
    CBF['RAW'] = {}
    CBF['RAW']['OTHERPRIORS']=np.expand_dims(OPR,axis=1)
    CBF['RAW']['OTHERPRIORSUNC']=np.expand_dims(OPRU,axis=1);
    CBF['RAW']['info']='Raw inputs/outputs as stored in CBF binary structure';
    CBF['RAW']['details']='For completeness & development purpose only; When re-writing CBF to file, these are over-written by CBF.OBS, etc.';

    
    
    
    
    return CBF
    #disp(sprintf('CHECK: .cbf file "%s" successfully read into matlab.',filename))   


def read_cbr_file(filename,INFO):

    """
    Read binary cbr file with parameters in it
    
    Consider expanding this function to read multiple files of the same
    run.
    
    Adapted from readbinarymat.m and CARDAMOM_READ_BINARY_FILEFORMAT.m from
    Anthony Bloom
    
    readbinarymat.m is set up to handle binary files with more than two 
    dimensions.
    
    Not sure what 'latterhalf' number is for
    
    """
    
    
    # Defaul setting
    if not len(INFO):
        # Number of parameters
        INFO = {'nopars':32,
                'latterhalf':0}
    
    with open(filename, 'rb') as fid:
        BD = np.fromfile(fid, np.float64)

    N = int(BD.shape[0]/INFO['nopars'])
    
    dims = [N,INFO['nopars']]
    
    PARS = BD.reshape(dims)
    
    return PARS


def write_cbr_file(PARS, outfilename):
    """
    INPUTS: 
        - PARS: CARAMOM parameter array Ensemble X Parameter
        - "outfilename": name of file 
        Writes PARS array to CARDAMOM binary file (.cbr) format
        
    """
    dims = PARS.shape
    BINARY_DATA = PARS.reshape(dims[0]*dims[1])

    
    # Write data to binary
    with open(outfilename, 'wb') as fid:
        BINARY_DATA.tofile(fid)



def CARDAMOM_READ_BINARY_FILEFORMAT(filename,cbrfrac=0.5,INFO=[]):
    """
    Take either a .cbf file, or a .cbr file, or a list of .cbr files
    """
    
    if type(filename) is str:
        
        if filename.split(".")[-1] == "cbf":
            VARRETURN = read_cbf_file(filename)
        
        elif filename.split(".")[-1] == "cbr":
            PARS = read_cbr_file(filename,INFO)
            VARRETURN = PARS[int((PARS.shape[0])*cbrfrac):]
            
            
    elif type(filename) is list:
        
        filename.sort()
        tmplist = []
        allcbrlist = []
        for fn in filename:
            
            if fn.split(".")[-1] == "cbr":
                PARS = read_cbr_file(fn,INFO)
                tmplist.append(PARS[int((PARS.shape[0])*cbrfrac):])
                allcbrlist.append(True)
            else:
                allcbrlist.append(False)
        
        VARRETURN = np.concatenate(tmplist,axis=0)
        
        if not all(allcbrlist):
            print('not all files where .cbr files')
            
    return VARRETURN
               


def compile_cbf_obs_fields(CBF):
    '''
    %Time-resolved observations
    '''
    
    if 'OBS' in CBF.keys():
    
        if type(CBF['OBS']) == np.ndarray:
            CBFOBS=CBF['OBS'];
            warnings.warn('CARDAMOM_WRITE_BINARY_FILEFORMAT: CBF.OBS matrix used to define observations (will soon be obsolete)', DeprecationWarning)
        
        else:
        
            obsnames=['GPP','LAI','NBE','ABGB','ET','EWT','BAND1','BAND2','BAND3','BAND4','SOM']
            fnames=CBF['OBS'].keys()
            
            isobs = []
            for obn in obsnames:
                isobs.append(obn in fnames)
            
            for jjob, tfobs in enumerate(isobs):
                if tfobs == True:
                    lastob = jjob
            nobs = lastob+1     
            
            #Defining CBFOBS
            N=np.shape(CBF['MET'])[0]
            CBFOBS=np.ones((N,np.max([nobs,3])))*-9999.0
            
            for n in np.arange(len(isobs))[isobs]:
                #Empty fields are tolerated (consequently no writing to file)
                if not not list(CBF['OBS'][obsnames[n]]) and np.shape(CBF['OBS'][obsnames[n]])[0]== N:
                    CBFOBS[:,n]=CBF['OBS'][obsnames[n]]
                if np.shape(CBF['OBS'][obsnames[n]])[0]!= N and not not list(CBF['OBS'][obsnames[n]]):
                    warnings.warn('CBF[OBS] field dimensions not compatible with CBF.MET, writing -9999')

    
    else:
        #FIlling in OBS with -9999 values;  
        #Eventually step can be made obsolete when C code can accept
        #(1) noobs, and (2) obsid
         
        CBFOBS=np.ones((np.shape(CBF['MET'])[0],3))*-9999
    
    return CBFOBS




def write_obs_uncertainty_fields(CBF,SD,OPRU):
    '''
    %Time-resolved observation uncertainties
    UPDATE IN PROGRESS
    
    %From CARDAMOM_READ_BINARY_DATA.c 
    % DATA->nee_annual_unc=statdat[13];
    % DATA->et_annual_unc=statdat[14];
    % DATA->nee_obs_unc=statdat[15];if (statdat[15]<0){DATA->nee_obs_unc=0.5;}
    % DATA->et_obs_unc=statdat[16];if (statdat[16]<0){DATA->et_obs_unc=2;}
    % DATA->ewt_annual_unc=statdat[17];
    % DATA->ewt_obs_unc=statdat[18];if (statdat[18]<0){DATA->ewt_obs_unc=50;}
    % DATA->gpp_annual_unc=statdat[19];
    % DATA->gpp_obs_unc=statdat[20];if (statdat[18]<0){DATA->gpp_obs_unc=2;}
    '''



    #NBE
    SD[13]=CBF['OBSUNC']['NBE']['annual_unc']
    SD[15]=CBF['OBSUNC']['NBE']['seasonal_unc']
    #CBF.OBSUNC.NBE.info='Single point (default = 0.5) and annual (annual_unc) NBE uncertainty [gC/m2/d]';
    
    #ET 
    SD[14]=CBF['OBSUNC']['ET']['annual_unc']
    SD[16]=CBF['OBSUNC']['ET']['unc']
    #CBF.OBSUNC.ET.info='Single point (default = 2) and annual (annual_unc) ET uncertainty [mm/d]';
    
    #EWT
    SD[17]=CBF['OBSUNC']['EWT']['annual_unc']
    SD[18]=CBF['OBSUNC']['EWT']['unc']
    #CBF.OBSUNC.EWT.info=sprintf('Single point (default = 50) and annual (N/A yet) EWT uncertainty [mm]');
    
    #GPP 
    SD[19]=CBF['OBSUNC']['GPP']['annual_unc']
    SD[20]=CBF['OBSUNC']['GPP']['unc']
    #CBF.OBSUNC.GPP.info='Single point (default = 2) and annual (annual_unc) GPP uncertainty [gC/m2/d]';
    #gpp abs: option on treating GPP constraint as absolute or relative
    SD[7]=CBF['OBSUNC']['GPP']['gppabs']
    #CBF.OBSUNC.GPP.gppabs_info='Set to "1" for GPP data, set to "0" for SIF data"';
    
    #Threshold values (for positive-definite values)
    SD[21]=CBF['OBSUNC']['ET']['obs_unc_threshold'];
    SD[22]=CBF['OBSUNC']['GPP']['obs_unc_threshold'];
    
    
    
    #Time-resolved biomass uncertaiknty
    OPRU[1]=CBF['OBSUNC']['ABGB']['unc']
    #CBF.OBSUNC.ABGB.info=sprintf('Uncertainty on time-resolved biomass CBF.OBS.ABGB');
    #Time-resolved SOM uncertainty
    OPRU[7]=CBF['OBSUNC']['SOM']['unc']
    #CBF.OBSUNC.SOM.info=sprintf('Single point uncertainty on time-resolved SOM, CBF.OBS.SOM');
    
    
    #Note: all "otherpriors" and "otherpriorsunc" will be passed as follows
    if type(CBF['OBS']) == dict and 'ABGB' in CBF['OBS'].keys() and not not list(CBF['OBS']['ABGB']) and CBF['OBSUNC']['ABGB']['unc']==-9999:
        print('Need to prescribe CBF.AGBGunc field')
    if type(CBF['OBS']) == dict and 'SOM' in CBF['OBS'].keys() and not not list(CBF['OBS']['SOM']) and CBF['OBSUNC']['SOM']['unc']==-9999:
        print('Need to prescribe CBF.SOMunc field')

    return [SD,OPRU]



 
def write_other_obs_constraints(CBF,OPR,OPRU):
    '''
    %other (time-invariant) constraints
    '''

    #Mean biomass
    OPR[0]=CBF['OTHER_OBS']['MBiomass']['mean']
    OPRU[0]=CBF['OTHER_OBS']['MBiomass']['unc']
    
    #Mean fire emissions
    OPR[2]=CBF['OTHER_OBS']['MFire']['mean']
    OPRU[2]=CBF['OTHER_OBS']['MFire']['unc']
    
    #Mean LAI
    OPR[4]=CBF['OTHER_OBS']['MLAI']['mean']
    OPRU[4]=CBF['OTHER_OBS']['MLAI']['unc']
    
    #Mean GPP
    OPR[5]=CBF['OTHER_OBS']['MGPP']['mean']
    OPRU[5]=CBF['OTHER_OBS']['MGPP']['unc']

    return [OPR,OPRU]


    
def CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF,outfilename):
    
    """
    INPUTS: 
        - CBF: CARAMOM binary structure
        - "filename": name of file 
        Writes CBF structure to CARDAMOM binary file (.cbf) format
        See CARDAMOM_READ_BINARY_FILEFORMAT.m for obtaining a "CBF" structure
        template
        
        Adapted in 2018/04/17 from matlab script by A.A. Bloom 2018/02/11
        %Last updated from A.A. Bloom 2019/07/24
    """
    
    fill_value = -9999
    SD = np.ones((100,1))*fill_value
    
    PR = CBF['PARPRIORS']
    PRU = CBF['PARPRIORUNC']
    OPR=CBF['RAW']['OTHERPRIORS']
    OPRU=CBF['RAW']['OTHERPRIORSUNC']
    

    #Number of timesteps
    nodays=np.shape(CBF['MET'])[0]
    #Number of met fields
    nomet=np.shape(CBF['MET'])[1]
    #
    CBF['OBS']=compile_cbf_obs_fields(CBF)
    #Number of time-resolved observations
    noobs=np.shape(CBF['OBS'])[1]
    
    CBF['nomet']=nomet;
    CBF['noobs']=noobs;
    CBF['nodays']=nodays;
    
    if CBF['nomet']!=nomet:
        warnings.warn('CBF.nomet updated for consistency with CBF.MET dimensions')
        CBF['nomet']=nomet
    if CBF['noobs']!=noobs:
        warnings.warn('CBF.noobs updated for consistency with CBF.OBS dimensions')
        CBF['noobs']=noobs
    if CBF['nodays']!=nodays:
        warnings.warn('CBF.nodays updated for consistency with CBF.MET and CBF.OBS dimensions')
        CBF['nodays']=nodays
    
    # Must update if new information is available
    SD[0:7] = np.expand_dims(np.array([CBF['ID'],
                                      CBF['LAT'],
                                      CBF['nodays'],
                                      CBF['nomet'],
                                      CBF['noobs'],
                                      CBF['EDC'],
                                      CBF['EDCDIAG']]),axis=1)
    
    
    
    # EDCDIAG - set "on" by default
    # if isfield(MD,'EDCDIAG')==0;SD(7)=1;end

    #Write time-invariant terms
    SD,OPRU=write_obs_uncertainty_fields(CBF,SD,OPRU)

    #Write other obs constraints
    OPR,OPRU=write_other_obs_constraints(CBF,OPR,OPRU)
    
    
    # Prescribe reference met
    if 'mmet' in CBF.keys() and CBF['mmet'].shape[0] == CBF['nomet']:

        MTEMPDATA = np.expand_dims(
                np.concatenate((CBF['mmet'],np.ones((CBF['noobs']))*fill_value)),
                axis=0)
            
        # Order = timestep1 (met, obs), timestep2 (met, obs)....
        TSDATA = np.concatenate((CBF['MET'],CBF['OBS']),axis=1)
        TEMPDATA = np.expand_dims(np.concatenate((TSDATA,MTEMPDATA)).flatten(),axis=1)
            
    else:
        MTEMPDATA = np.array([],ndmin=2)
        
        # Order = timestep1 (met, obs), timestep2 (met, obs)....
        TSDATA = np.concatenate((CBF['MET'],CBF['OBS']),axis=1)
        TEMPDATA = np.expand_dims(TSDATA.flatten(),axis=1)
    
        
    if TEMPDATA.shape[0] != (CBF['nomet']+CBF['noobs'])*CBF['nodays']+MTEMPDATA.shape[1]:
        print('Error!! incorrectly set OBS or MET data vector!!')
    else:
        BINARY_DATA = np.concatenate((SD,PR,PRU,OPR,OPRU,TEMPDATA))
    
    # Write data to binary
    with open(outfilename, 'wb') as fid:
        BINARY_DATA.tofile(fid)

def CARDAMOM_READ_OUTPUT(cbffile,cbrfile,fluxfile,poolfile,probfile=[],baddata = -9999,INFO=[]):
    
    
    # Defaults for ID = 803 - DALEC_FIREBUCKET3
    CBR = {}
    
    nofluxes = 30
    nopools = 7
    CBF = read_cbf_file(cbffile)
    nopars = len(readcardamommodel.get_parnames(int(CBF['ID']),output='shortnames'))
    model_info = {'nopars':nopars,
                  'latterhalf':0}
    
    PARS = read_cbr_file(cbrfile,INFO = model_info)
    nodays = CBF['nodays']

    FLUXES = readbinarymat(fluxfile,[nodays,nofluxes])
    if FLUXES is None:
        return None
    
    POOLS = readbinarymat(poolfile,[nodays+1,nopools])
    if len(probfile) > 0:
        PROBS = readbinarymat(probfile,[1])
    
    FLUXES[FLUXES == baddata] = np.nan
    POOLS[POOLS == baddata] = np.nan
    if len(probfile) > 0:
        PROBS[PROBS == baddata] = np.nan
    
    CBR['FLUXES'] = FLUXES
    CBR['POOLS'] = POOLS[:,1:,:]
    if len(probfile) > 0:
        CBR['PROBS'] = PROBS
    
    # GPP
    CBR['GPP'] = FLUXES[:,:,0]
    
    CBR['NEE'] = np.sum(FLUXES[:,:,[2,12,13]],axis=2) - CBR['GPP']
    
    # NBE
    if CBF['ID'] > 1:
        CBR['NBE'] = CBR['NEE'] + CBR['FLUXES'][:,:,16]
    else:
        CBR['NBE'] = CBR['NEE']
    
    # LAI
    if CBF['ID'] == 821:
        # Shift location of LMCA for calculation in 821
        CBR['LAI']=CBR['POOLS'][:,:,1]/np.expand_dims(PARS[:,15],1)
    else:
        CBR['LAI']=CBR['POOLS'][:,:,1]/np.expand_dims(PARS[:,16],1)
    
    #Water Stress [something to add]
#        if CBR['POOLS'].shape[-1]>6:
#            if CBF['ID']<= 8 or any([CBF['ID']==s for s in [801,802,803]]):
#                CBR['H2OSTRESS'] = 

    return CBR 

def CARDAMOM_READ_OUTPUT_FILEAWARE(cbffile,cbrfiles,outputfiles,cutfrac = 0.5,INFO=[],chaininfoout=None):
    
    # Output names present
    namedvars = list(set([outf.split("/")[-1].split("_")[0] for outf in outputfiles]))
    
    headdir = "/".join(outputfiles[0].split("/")[:-2]+[""])
    headdir_cbr = "/".join(cbrfiles[0].split("/")[:-2]+[""])
    cbrfilelist = cbrfiles
    fluxfilelist = [op for op in outputfiles if 'fluxfile' in op]
    poolfilelist = [op for op in outputfiles if 'poolfile' in op]
    
    cbrname = "_".join(cbrfiles[0].split("/")[-1].split(".")[0].split("_")[:-1])
    fluxname = "_".join(fluxfilelist[0].split("/")[-1].split(".")[0].split("_")[:-1])
    poolname = "_".join(poolfilelist[0].split("/")[-1].split(".")[0].split("_")[:-1])
    
    if 'probfile' in namedvars:
        probfilelist = [op for op in outputfiles if 'probfile' in op]
        probname = "_".join(probfilelist[0].split("/")[-1].split(".")[0].split("_")[:-1])
    
    
    #
    ## Check for which runs they each have
    count = 0
    runnums = []
    
    
    if 'probfile' in namedvars:
        runlist = [cbrfilelist,fluxfilelist,poolfilelist,probfilelist]
    else:
        runlist = [cbrfilelist,fluxfilelist,poolfilelist]
    
    while count < len(runlist)-1:
        if count == 0:
            runnums = [s.split("/")[-1].split(".")[0].split("_")[-1] for s in runlist[count]]
            
        runnums = list(set(runnums) & 
                       set([s.split("/")[-1].split(".")[0].split("_")[-1] for s in runlist[count+1]]))
        count = count + 1
    
    runnums = list(np.sort(runnums)) # Keep in the same order as other files
    tmp = {}
    Nlist = []
    for idx,rn in enumerate(runnums):
        
        cbrfile = headdir_cbr + 'cbr/'+cbrname+'_'+rn+'.cbr'
        fluxfile = headdir + 'output/'+fluxname+'_'+rn+'.bin'
        poolfile = headdir + 'output/'+poolname+'_'+rn+'.bin'
        
        if 'probfile' in namedvars:
            probfile = headdir + 'output/'+probname+'_'+rn+'.bin'
    
            CBRtmp = CARDAMOM_READ_OUTPUT(cbffile,cbrfile,fluxfile,poolfile,probfile,INFO)
            
        else:
            CBRtmp = CARDAMOM_READ_OUTPUT(cbffile,cbrfile,fluxfile,poolfile,INFO)
        
        if CBRtmp is not None:
            N = int(cutfrac*CBRtmp['FLUXES'].shape[0])
            Nlist.append(CBRtmp['FLUXES'].shape[0] - N)
        else:
            continue
        
        for ky in CBRtmp:
            if idx == 0:
                tmp[ky] = []
            
            if len(CBRtmp[ky].shape) == 1:
                tmp[ky].append(CBRtmp[ky][N:])
            elif len(CBRtmp[ky].shape) == 2:    
                tmp[ky].append(CBRtmp[ky][N:,:])
            
            elif len(CBRtmp[ky].shape) == 3:
                tmp[ky].append(CBRtmp[ky][N:,:,:])
            
    CBRall = {}
    for ky in tmp:
        CBRall[ky] = np.concatenate(tmp[ky],axis=0)
    
    if chaininfoout == 'list':
        return CBRall, Nlist 
    else:
        return CBRall
        

# %%

if __name__ == '__main__':

    # .cbf files
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/erai_obs_2010_2015/cbf/'
    fnbinary = 'GCRUN_NOV17_20_1623.cbf'
    fnbinaryout = 'GCRUN_NOV17_20_1623_testtest.cbf'
    
    # .cbr files
    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/erai_obs_2010_2015/cbr/'
    fncbrbinary = 'GCRUN_NOV17_20_1623_1.cbr'
    
    # .cbr file list
    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/erai_obs_2010_2015/cbr/'
    fncbrbinarylist = 'GCRUN_NOV17_20_1623*.cbr'
    
    #
    # Read a .cbf file into a CBF dictionary
#    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p12_uncmatch_lens001/cesmlens/cbf/'
#    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p12_uncmatch_full/cruncep/cbf/'
#    fnbinary = 'cruncep001GCRUN_NOV127_20_smoothnbe_sifnojflwlai_1997_2016_nbeseaunc100_100_nbeannunc100_2_1623.cbf'
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p16_long/cruncep_osse1/cbf/'
    fnbinary = 'cruncep001GCRUN_ossemedian_1526.cbf'
#    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p17_long_test3/cruncep/cbf/'
#    fnbinary = 'cru004GCR003_1997_2016_nbe2002_2042.cbf'
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p18_long/cruncep/cbf/'
    fnbinary = 'cru004GCR005_1920_2016_nbe2002_2042.cbf'

    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p18_long/cesmlensbc001/cbf/'
    fnbinary = 'lbc001GCR005_1920_2016_nbe2002_1623.cbf'
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p18_long/lnsbcco2no001/cbf/'
    fnbinary = 'bcco2n001GCR005_1920_2016_nbe2002_2042.cbf'
    
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/p20_long/cruncep/cbf/'
    fnbinary = 'cru004GCR006_1920_2015_nbe2002_2042.cbf'
    
    #
    ## Test file
    dirbinary = MACHDIR + '/Google Drive/DATA/DALEC/harvard_test/lai_bias/cbf/'
    fnbinary = 'e1a_USHa1_coplai_div100_lwbias100_0000.cbf'
    
    inputfilename = dirbinary+fnbinary
    CBF1 = read_cbf_file(inputfilename)
    

    #
    # Write a binary .cbf file from a given CBF dictionary
    fnbinaryout = 'e1a_USHa1_id821_test_coplai_div100_lwbias100_0000.cbf' #'e1b_USHa1_821_0001.cbf'
    outfilename = dirbinary+fnbinaryout
    CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF_new,outfilename)
    
    
    #
    # Reading a .cbr file
#    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/p12_mbptest_819/cruncep/cbr/'
#    fncbrbinary = 'cruncep001GCRUN_NOV127_20_smoothnbe_sifnojflwlai_1997_2016_nbeseaunc100_100_nbeannunc100_2_3313_1.cbr'
    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/p18_long/cruncep/cbr/'
    fncbrbinary = 'cru004GCR005_1920_2016_nbe2002_2042_1.cbr'
    
    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/harvard_test/lai_bias/cbr/'
    fncbrbinary = 'e1a_USHa1_id821_test_coplai_div100_lwbias100_0000.cbr'
#    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/p17_long_test3/cruncep/cbr/'
#    fncbrbinary = 'cru004GCR003_1997_2016_nbe2002_3260_1.cbr'
#    dircbrbinary = MACHDIR + '/Google Drive/DATA/DALEC/p14_uncmatch/cruncep/cbr/'
#    fncbrbinary = 'cru004GCR003_1920_2016_nbe2002_3436_1.cbr'
    
    
    inputfilename_cbr = dircbrbinary + fncbrbinary
    
    PARS = read_cbr_file(inputfilename_cbr,INFO = {'nopars':36,'latterhalf':0})
    
    #
    # Reading a list of .cbr files
    inputfilename_cbrlist = glob.glob(dircbrbinary + fncbrbinarylist)
    PARS = CARDAMOM_READ_BINARY_FILEFORMAT(inputfilename_cbrlist)
    
    cbrtestfilename = dircbrbinary + 'testtesttest.cbr'
    write_cbr_file(PARS, cbrtestfilename)
    
    
    
    
    #
    # Test write out of CBF, PARS to matlab .mat for matlab DALEC
    testfile = MACHDIR + '/Google Drive/DATA/DALEC/test/cbr/testcbrreadwrite'
    sio.savemat(testfile+'.mat',{'CBF':CBF,'PARS':PARS})
    

 
    #%% Open output files
    testdir = MACHDIR + '/Google Drive/DATA/DALEC/p20_long/cruncep/output/'
    testcbfdir = MACHDIR + '/Google Drive/DATA/DALEC/p20_long/cruncep/cbf/'
    testcbrdir = MACHDIR + '/Google Drive/DATA/DALEC/p20_long/cruncep/cbr/'
    
    cbfname = 'cru004GCR006_1920_2015_nbe2002_2042.cbf'
    cbffile = testcbfdir + cbfname
    
    cbrname = 'cru004GCR006_1920_2015_nbe2002_2042_3.cbr'
    cbrfile = testcbrdir + cbrname
    
    fluxname = 'fluxfile_cru004GCR006_1920_2015_nbe2002_2042_3.bin'
    fluxfile = testdir + fluxname
    poolname = 'poolfile_cru004GCR006_1920_2015_nbe2002_2042_3.bin'
    poolfile = testdir + poolname
    edcdname = 'edcdfile_cru004GCR006_1920_2015_nbe2002_2042_3.bin'
    edcdfile = testdir + edcdname
    probname = 'probfile_cru004GCR006_1920_2015_nbe2002_2042_3.bin'
    probfile = testdir + probname
    
    probfilemat = readbinarymat(probfile,[1])
    
    edcdfilemat = readbinarymat(edcdfile,[250,100])
    
#    with open(edcdfile, 'rb') as fid:
#        BD = np.fromfile(fid, np.int32)
    
    CBR = CARDAMOM_READ_OUTPUT(cbffile,cbrfile,fluxfile,poolfile)
    CBRall = CARDAMOM_READ_OUTPUT_FILEAWARE(cbffile)
    
    CBF = read_cbf_file(cbffile)
    for cbf in CBF:

        if type(CBF[cbf]) == list:
            nshape = len(CBF[cbf])
        elif type(CBF[cbf]) == np.ndarray:
            nshape = CBF[cbf].shape
        else:
            nshape = 1

        print(cbf+' , '+str(nshape))        
            
    for cbf in CBF:
        
        if type(CBF[cbf]) != list and type(CBF[cbf]) != np.ndarray:
            print(cbf +','+str(CBF[cbf]))   
            
    #%% Testing the new .cbr synthetic files
    testcbrdir = MACHDIR + '/Google Drive/DATA/DALEC/p2_uncmatch/erainterim/cbr/'
    cbrname = 'CBFerai_2001-2016_obsmod_icabg_fmodis_clstr_11_divide_median_rndm_2042_0.cbr'
    cbrfile = testcbrdir + cbrname
    PARS1 = read_cbr_file(cbrfile,[])

    PARS1_mean = PARS1.mean(axis=0)
    PARS_mean = PARS.mean(axis=0)            
    
    #%% Check newest .cbf files
    
    testcbfdir = MACHDIR + '/Google Drive/DATA/DALEC/p5_uncmatch_v1/erainterim/cbf/'
    
    
    cbfname = 'CBFerai_2001-2016_obsab_icabg_fmodis_uncmatch_4359.cbf'
    cbffile = testcbfdir + cbfname
    CBF = read_cbf_file(cbffile)       
    
    
   
