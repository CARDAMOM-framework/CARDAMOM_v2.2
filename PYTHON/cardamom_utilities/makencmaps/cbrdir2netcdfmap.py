#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:02 2020

@author: gquetin


Projects: J5, J6
Convert a dictionary of cbr bin files from CARDAMOM to a map/netcdf of parameters with names etc.

Additional function of creating a directory of cbr bin files from the netcdf file

"""


#import numpy as np
import xarray as xr

import os
import glob
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 

import loadcmip5todalecmet as ldcmip5
import readcardamommodel
import analycbr

imp.reload(ldcmip5)



# %%

if __name__ == '__main__':
    
    plotsave = False
    datasave = True
    mapplot = False
    makeplots = False
    
    # Set up folders based on the home folder
    datadir, dataoutdir = ldcmip5.folderdeff(MACHDIR)
    
    exper = 'p20_long' #'p16_long' #'p6_uncmatch_v1_smthnbe'
    model = 'cruncep' #'erainterim'
    cutfrac = .5 # may have this labeled wrong, says how much to keep rather than throw away
    filesize = 144000 # 144000
    
    dircbr = MACHDIR + '/Google Drive/DATA/DALEC/'+exper+'/'+model+'/cbr/'
    dircbr_original = MACHDIR + '/Google Drive/DATA/DALEC/erai_obs_2010_2015/cbr/'
    dircbrout = MACHDIR + '/Google Drive/DATA/DALEC/'+exper+'/'+model+'/cbr/'
    plotdir = MACHDIR + '/Google Drive/DATA/Analysis/DALEC_CMIP5/Models/erainterim/Plots/summary/'
    
    
#    dircbr = '/Users/gquetin/Google Drive/CARDAMOM_Share/SharedNotes_Experiments/ball_berry_model_tests/data/cbr/'
    
    dictfiles = {
            'GCRUN0':{'wildcard':dircbr_original+'*GCRUN*.cbr',
                      'filelist':[]},
            'GCRUN1':{'wildcard':dircbr+'*GCRUN*.cbr',
                      'filelist':[]},
            'GCRUN_fireclim':{'wildcard':dircbr+'*mclim_*.cbr',
                      'filelist':[]},
            'GCRUN_fireclimmod':{'wildcard':dircbr+'*mclimmod*.cbr',
                      'filelist':[]},
            'p3_uncmatch_limit_20':{'wildcard':dircbr + '*2001-2016_obsab_icabg_fmodis_uncmatch*.cbr'},
            'p2_uncmatch':{'wildcard':dircbr+'*2001-2016_obsmod_icabg_fmodis_uncmatch*.cbr'},
            'p6_smthnbe_raw':{'wildcard':dircbr+'*2001-2016_obsmod_sif20102015lai20102015nbe20102013mlaion_icabg_fmodis_uncmatch*.cbr'},
            'p8_uncmatch':{'wildcard':dircbr+'cruncep*.cbr'},
            'p9_uncmatch_sea10_ann20':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeunc100_10_*.cbr'},
            'p9_uncmatch_sea50_ann20':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeunc100_50_*.cbr'},
            'p9_uncmatch_sea100_ann20':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeunc100_100_*.cbr'},
            'p9_uncmatch_sea150_ann20':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeunc100_150_*.cbr'},
            'p9_uncmatch_sea200_ann20':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeunc100_200_*.cbr'},
            'p10_uncmatch_sea50_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeseaunc100_50_nbeannunc100_2_*.cbr'},
            'p10_uncmatch_sea100_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeseaunc100_100_nbeannunc100_2_*.cbr'},
            'p10_uncmatch_sea200_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeannunc100_2_*.cbr'},
            'p11_uncmatch_sea100_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_lwlaisif_1997_2016_nbeseaunc100_100_nbeannunc100_2_*.cbr'},
            'p12_uncmatch_sea100_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_sifnojflwlai_1997_2016_nbeseaunc100_100_nbeannunc100_2_*.cbr'},
            'p13_uncmatch_sea200_ann2':{'wildcard':dircbr+'cruncep001GCRUN_NOV127_20_smoothnbe_sifnojflwlai_1997_2016_nbeseaunc100_200_nbeannunc100_2_*.cbr'},
            'p16_long':{'wildcard':dircbr+'cru004GCR005_1920_2016_nbe2002_*.cbr'},
            'p16_long_osse':{'wildcard':dircbr + 'cru004GCR005_1920_2016_nbe2002_ossemedian_*.cbr'},
            'p16_long_osse4':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p16_long_osse1':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p16_long_osse2':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p16_long_osse3':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p16_long_osse3b':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p16_long_osse5':{'wildcard':dircbr + 'cruncep001GCRUN_ossemedian_*.cbr'},
            'p12_carolinebbrun':{'wildcard':dircbr + '821_cruncep001GCRUN_NOV127_20_smoothnbe_sifnojflwlai_1997_2016_nbeseaunc100_100_nbeannunc100_2_*.cbr'},
            'p17_long':{'wildcard':dircbr + 'cru004GCR005_1920_2016_nbe2002_*.cbr'},
            'p18_long_temp':{'wildcard':dircbr + 'cru004GCR005_1920_2016_nbe2002_*.cbr'},
            'p18_long_temp2':{'wildcard':dircbr + 'cru004GCR005_1920_2016_nbe2002_*.cbr'},
            'p19_shorttest':{'wildcard':dircbr + 'cru004GCR005_1999_2015_nbe2002_*.cbr'},
            'p20_test':{'wildcard':dircbr + 'cru004GCR005_1920_2016_nbe2002_*.cbr'},
            'p20_shorttest':{'wildcard':dircbr + 'cru004GCR006_1999_2015_nbe2002_*.cbr'},
            'p20_long':{'wildcard':dircbr + 'cru004GCR006_1920_2015_nbe2002_*.cbr'},
            'p21_luctest':{'wildcard':dircbr + 'cru004GCR006_1920_2015_nbe2002_luc_*.cbr'}
            }
    
    modelid = 821 #809 
    filepicklist = ['p20_long']

    
    
    pnames_org_short = readcardamommodel.get_parnames(modelid,output='shortnames')
            
    INFO = {'nopars':len(pnames_org_short),
            'latterhalf':0}
#    dircbr = MACHDIR + '/Google Drive/DATA/DALEC/multipt_test/HadGEM2-ES/cbr/'
#    plotdir = MACHDIR + '/Google Drive/DATA/Analysis/DALEC_CMIP5/multipt_test/HadGEM2-ES/Plots/'
#    
#    parrunlist = ['GCRUN','2005mclimmod','2005mclim']
#    wildcard = [dircbr1+'*GCRUN*.cbr',dircbr + '*2005-2015*mclimmodel*.cbr',dircbr + '*2005-2015*mclim_*.cbr']
    lparnames = pnames_org_short #list_parnames(modelid=modelid)
    
    
    
    
    
    dsdict = {}
    for dfile in filepicklist:
        
        dictfiles[dfile]['filelist'] = glob.glob(dictfiles[dfile]['wildcard'])
        #
        # Put in a list of points here to limit
        fl = dictfiles[dfile]['filelist']
        
        
        dircbrhere = "/".join(fl[0].split("/")[:-1]+[''])
        
        # Save out to netCDF
        savename = dircbrhere + "_".join(['PARS',fl[0].split("/")[-1][:-11],'all'])
        if (os.path.exists(savename+'.nc')):
            # Load from file
            dsdict[dfile] = xr.open_dataset(savename+'.nc')
            
        else:
            dsdict[dfile] = analycbr.cbr_to_netcdf(fl,cutfrac=cutfrac,INFO=INFO,filesize=filesize,modelid=modelid)
            
            if datasave is True:
                dsdict[dfile].to_netcdf(savename+'.nc')  
        
    
    if 'GCRUN1' in filepicklist:
        dsdict['GCRUN1_match'] = dsdict['GCRUN1'].sel(lat=slice(dsdict['GCRUN0'].lat.values.min(),
                                                dsdict['GCRUN0'].lat.values.max()),
                                        lon=slice(dsdict['GCRUN0'].lon.values.min(),
                                                dsdict['GCRUN0'].lon.values.max()))

    
    #%
    dsdict_med = {}
    dsdict_std = {}
    
    for dd in dsdict:
        
        dsdict_med[dd] = dsdict[dd].median(dim='ensemble').copy()
        dsdict_std[dd] = dsdict[dd].std(dim='ensemble').copy()
    
    #%% Write a binary .cbr from a netcdf file (single point) and then (multiple point)
    
    testwritebinaryflag = False
    
    if testwritebinaryflag == True:
        
        ftag = dictfiles[dfile]['wildcard'].split("/")[-1].split('*')[0][:-1]
        data_set_xy = dsdict[filepicklist[0]]
        
        yesno_list = analycbr.netcdfmap_to_cbrs(data_set_xy,nametag=dircbrout+ftag,modelid = modelid)  
        
#    ############# ALL PLOT SCRIPTS FROM HERE ###############
#    
#    
#    
#    
#    #%% Plot maps of difference of JPL out
#    difmapplot = False
#    if difmapplot == True:
#        parnamearray = np.stack(np.split(np.array(parnames(modelid=modelid,output='shortnames')),16))
#        vararray = np.stack(np.split(np.arange(1,33),16))
#        if mapplot == True:
#            
#            import matplotlib.pyplot as plt
#            import mpl_toolkits.basemap as bm #This has a problem with 'pyproj'
#            from mpl_toolkits.basemap import Basemap
#            
#            parnamearray = np.stack(np.split(np.array(parnames(modelid=modelid,output='longnames')),8))
#            vararray = np.stack(np.split(np.arange(1,33),8))
#            
#            diffcompare = ['GCRUN0','GCRUN1_match']
#            
#            #
#            # Mean maps
#            for varc,pnames in zip(vararray,parnamearray):
#                fig, axarr = plt.subplots(2,2,figsize=(12,9))
#                axarrf = axarr.flatten()
#                
#                pnlist = []
#                for vc,ax,pn in zip(varc,axarrf,pnames[:4]):
#            
#            
#                    #plot data
#                    dsref = dsdict[diffcompare[0]]
#                    var = pn
#                    
#                    pnlist.append(str(vc))
#                    
#                    if var != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*':
#                        pltdiff = (dsdict_med[diffcompare[0]][var].values - dsdict_med[diffcompare[1]][var].values)
#                        pltstd = (dsdict_std[diffcompare[0]][var].values + dsdict_std[diffcompare[0]][var].values)/2
#                        
#                        pltdata = pltdiff/pltstd #/ds_obs0med[var].values
#                        drange = 1.5*np.nanstd(pltdata)
#                        #pltdata = ds[var].sel(time='2010').values - ds[var].sel(time='2100').values
#                        
#                        lons, lats = np.meshgrid(dsref.lon.values,dsref.lat.values)
#                        
#                        #create masked array where nan=mask. pcolormesh does not like NaNs.
#                        data_mask = np.ma.masked_where(np.isnan(pltdata.squeeze()),pltdata.squeeze())
#                        
#                        
#            
#                        #ax1 = fig.add_axes([0.05,0.05,0.9,0.9])
#                        #add title
#                        title_name = '(AB-GQ)/AB: ' + var
#                        #title_name = var +' 2010 - Wood C 2100'
#                        ax.set_title(title_name,fontsize=14)
#                        
#                        #cmap = plt.cm.BrBG
#                        cmap = plt.cm.PuOr
#                        #cmap = plt.cm.Reds
#                        #cmap = plt.cm.YlGnBu
#                        
#                        #create instance of Base map centered at 180
#                        centerlon = 0;
#                        m = Basemap(projection='moll',lon_0=centerlon,lat_0 = 0,resolution='c',ax=ax) #'l' is for low res.,'c' is for corse
#                        
#                        
#                        im1 = m.pcolormesh(lons,lats,data_mask,cmap=cmap,latlon=True,vmin=-drange,vmax=drange)
#                        
#                        #add lon and lat lines
#                        m.drawparallels(np.arange(-90.,90.,30.),labels=[1,1,0,0]) #labels = [left,right,top,bottom]
#                        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0])
#                        
#                        #draw map boundary
#                        m.drawcoastlines(linewidth=1.25)
#                        m.drawmapboundary(fill_color='1') #make map background white.To make it black, 
#                        
#                        #draw country borders
#                        m.drawcountries()
#                        
#                        # add colorbar
#                        cb = m.colorbar(im1,"bottom", size="5%", pad="15%")
#                    
#                plt.tight_layout()
#        
#            
#                if plotsave is True:    
#                    plotname = 'mapallpars_diffstd_'
#                                    
#                    savetag = "_".join(pnlist)
#                    plt.savefig(plotdir + plotname + savetag)
#                    #plt.savefig(plotdir + plotname + savetag +'.eps')
#                    plt.close("all")
#
#    #%% Plot the raw maps of parameters
#    if mapplot == True:
#
#        parnamearray = np.stack(np.split(np.array(parnames(modelid=modelid,output='longnames')),16))
#        vararray = np.stack(np.split(np.arange(1,33),16))
#        
#        diffcompare = filepicklist
#
#        #
#        # Mean maps
#        for varc,pnames in zip(vararray,parnamearray):
#            fig, axarr = plt.subplots(2,2,figsize=(12,9))
#            axarrf = axarr.flatten()
#            
#            
#                
#            pnlist = []
#            for vc,ax,pn in zip(varc,axarr,pnames):
#                pnlist.append(str(vc))
#                for axx,ky in zip(ax,diffcompare):
#            
#                    #plot data
#                    dsref = dsdict[diffcompare[0]]
#                    var = pn
#                    
#                    
#                    if var != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*':
#                        pltdata = dsdict[ky][var].mean(dim='ensemble').values
#                        dmean = np.nanmean(dsref[var].values)
#                        dstd = np.nanstd(dsref[var].values)
#                        drange = [0,dmean+2*dstd]
#                        
#                        
#                        lons, lats = np.meshgrid(dsref.lon.values,dsref.lat.values)
#                        
#                        #create masked array where nan=mask. pcolormesh does not like NaNs.
#                        data_mask = np.ma.masked_where(np.isnan(pltdata.squeeze()),pltdata.squeeze())
#                        
#                        
#            
#                        #ax1 = fig.add_axes([0.05,0.05,0.9,0.9])
#                        #add title
#                        title_name = ky + ': ' + var
#                        #title_name = var +' 2010 - Wood C 2100'
#                        axx.set_title(title_name,fontsize=14)
#                        
#                        #cmap = plt.cm.BrBG
#                        #cmap = plt.cm.PuOr
#                        #cmap = plt.cm.Reds
#                        cmap = plt.cm.YlGnBu
#                        
#                        #create instance of Base map centered at 180
#                        centerlon = 0;
#                        m = Basemap(projection='moll',lon_0=centerlon,lat_0 = 0,resolution='c',ax=axx) #'l' is for low res.,'c' is for corse
#                        
#                        
#                        im1 = m.pcolormesh(lons,lats,data_mask,cmap=cmap,latlon=True,vmin=drange[0],vmax=drange[1])
#                        
#                        #add lon and lat lines
#                        m.drawparallels(np.arange(-90.,90.,30.),labels=[1,1,0,0]) #labels = [left,right,top,bottom]
#                        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0])
#                        
#                        #draw map boundary
#                        m.drawcoastlines(linewidth=1.25)
#                        m.drawmapboundary(fill_color='1') #make map background white.To make it black, 
#                        
#                        #draw country borders
#                        m.drawcountries()
#                        
#                        # add colorbar
#                        cb = m.colorbar(im1,"bottom", size="5%", pad="15%")
#                
#            plt.tight_layout()
#    
#        
#            if plotsave is True:    
#                plotname = 'mappars_'
#                                
#                savetag = "_".join(diffcompare+pnlist)
#                plt.savefig(plotdir + plotname + savetag)
#                #plt.savefig(plotdir + plotname + savetag +'.eps')
#                plt.close("all")
#    
#
#    
#    #%% Compare histograms
#    
#    if makeplots == True:
#        
#        def drawbox_plot(data, offset,edge_color, fill_color):
#            pos = np.arange(data.shape[1])+offset 
#            bp = ax.boxplot(data, positions= pos, widths=0.3, patch_artist=True,
#                            showfliers=False,manage_xticks=False)
#            for element in ['boxes', 'whiskers', 'medians', 'caps']:
#                plt.setp(bp[element], color=edge_color)
#            for patch in bp['boxes']:
#                patch.set(facecolor=fill_color)
#        
#        from random import shuffle
#        
#        dschecklist = ['GCRUN0','GCRUN1_match']
#        refds = dsdict[dschecklist[1]].copy()
#        refmask = ~np.isnan(refds[list(refds.data_vars)[0]].mean(dim='ensemble'))
#        lons,lats = np.meshgrid(refmask['lon'].values, refmask['lat'].values)
#        llpts = list(zip(lats[refmask.values].flatten(),
#                        lons[refmask.values].flatten()))
#        ptlist = llpts.copy()
#        shuffle(ptlist)
#        ptlist = ptlist[0:2]
#        
#        for pt in ptlist:
#            fig, axarr = plt.subplots(8,4,figsize=(18,24))
#            axarrf = axarr.flatten()
#            
#            lat, lon = pt
#            for ax, pn in zip(axarrf,lparnames):
#        
#            
#                    
#                for fp in dschecklist:
#                    if pn != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*':
#                        
#                        
#                        pltdata = dsdict[fp][pn].sel(lat=lat,lon=lon).values
#        
#                    
#                        
#                        ax.hist(pltdata,alpha=.4,label=fp,density=True)
#                        axarrf[0].legend()
#        
#                    ax.set_title(pn)
#                    ax.set_ylabel('density')
#                    ax.set_xlabel(pn)
#                    
#                        
#                
#                plt.tight_layout()
#                
#                if plotsave is True:    
#                    plotname = 'histpars'
#                                    
#                    savetag = "_".join([plotname]+pnlist)
#                    plt.savefig(plotdir + plotname + savetag)
#                    #plt.savefig(plotdir + plotname + savetag +'.eps')
#                    plt.close("all")
#    
#    #%% Box plots of parameters at a point
#    
#    
#    
#        dschecklist = ['GCRUN0','GCRUN1_match']
#        refds = dsdict[dschecklist[1]].copy()
#        refmask = ~np.isnan(refds[list(refds.data_vars)[0]].mean(dim='ensemble'))
#        lons,lats = np.meshgrid(refmask['lon'].values, refmask['lat'].values)
#        llpts = list(zip(lats[refmask.values].flatten(),
#                        lons[refmask.values].flatten()))
#        ptlist = llpts.copy()
#        shuffle(ptlist)
#        ptlist = ptlist[0:2]
#    
#        parnamefix = [pn for pn in lparnames if pn != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*']
#        
#        for pt in ptlist:
#            lat, lon = pt
#            fig, ax = plt.subplots(figsize=(20,10))
#            colorpick = ["tomato","skyblue"]
#            shiftpick = [-0.2,0.2]
#            
#            for idx, (sp,cl,fp) in enumerate(zip(shiftpick,colorpick,dschecklist)):
#                boxds = dsdict[fp].sel(lat=lat,lon=lon)
#                
#                
#                datalist = []
#                for pn in parnamefix:
#                    datalist.append(boxds[pn].values)
#                    
#                dataarray = np.stack(datalist,axis=1)
#                if idx == 0:
#                    d0 = np.mean(dataarray,axis=0)
#            
#                drawbox_plot(dataarray/d0, sp, cl, "white")
#            
#            axes = plt.gca()
#    
#            axes.set_ylim([0,2])
#            
#            plt.xticks(np.arange(len(parnamefix)))
#            
#            plt.show()
#    
#        #%% Thumprints of the parameters (dot for median, shaded uncertainty)
#    
#        dschecklist = ['GCRUN0','GCRUN1_match']
#        refds = dsdict[dschecklist[1]].copy()
#        refmask = ~np.isnan(refds[list(refds.data_vars)[0]].mean(dim='ensemble'))
#        lons,lats = np.meshgrid(refmask['lon'].values, refmask['lat'].values)
#        llpts = list(zip(lats[refmask.values].flatten(),
#                        lons[refmask.values].flatten()))
#        ptlist = llpts.copy()
#        shuffle(ptlist)
#        ptlist = ptlist[0:64]
#        lats,lons = zip(*ptlist)
#        ptlist = np.array(ptlist)[np.argsort(list(lats))]
#    
#        parnamefix = [pn for pn in lparnames if pn != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*']
#        
#        # Pick standardization as median of all points of variable
#        normlist = []
#        for pn in parnamefix:
#            #normlist.append(np.nanmedian(dsdict[dschecklist[0]][pn].values.flatten()))
#            normlist.append(np.nanmedian(np.nanmedian(dsdict[dschecklist[0]][pn].values,axis=2).flatten()))
#            
#        d0 = np.array(normlist)  
#        
#        fig, axarr = plt.subplots(8,8,figsize=(40,20),sharex=True,sharey=True)
#        for ax,pt in zip(axarr.flatten(),ptlist):
#            lat, lon = pt
#            
#    
#            
#            for idx, fp in enumerate(dschecklist):
#                boxds = dsdict[fp].sel(lat=lat,lon=lon)
#                
#                dataarray = dataset_to_array(boxds,parnamefix)
#                
#                d75 = np.percentile(dataarray,75,axis=0)
#                d50 = np.percentile(dataarray,50,axis=0)
#                d25 = np.percentile(dataarray,25,axis=0)
#                
#                #if idx == 0:
#                    #d0 = np.median(dataarray,axis=0)
#            
#                ax.fill_between(np.arange(len(parnamefix)),d25/d0,d75/d0, alpha=0.4)
#                
#                
#                ax.plot(np.arange(len(parnamefix)),d50/d0,'o-',label=fp)
#            ax.set_title('lat: '+str(lat)+', lon: '+str(lon),fontsize=10)
#            #ax.set_xlabel('parameters',fontsize=8)
#            #ax.set_ylabel('xy std of median',fontsize=8)
#            ax.tick_params(labelsize = 8.0)
#            
#        plt.xticks(np.arange(len(parnamefix)))
#      
#        
#        axes = plt.gca()
#    
#        axes.set_ylim([0,8])
#        plt.tight_layout()
#        if plotsave is True:    
#            plotname = 'random64ptpar'
#                            
#            savetag = "_".join([plotname]+dschecklist)
#            plt.savefig(plotdir + savetag)
#            #plt.savefig(plotdir + plotname + savetag +'.eps')
#            plt.close("all")
#            
#            
#        #%% Pick regions
#        
#        
#        import summaryplots as sumplt
#        imp.reload(sumplt)
#        
#        regionlistall = sumplt.regionall_pickpolygon()
#        regionlisttropic = sumplt.regiontropic_pickpolygon()
#        
#        regionlist = regionlisttropic
#        datasettag = 'cbrregions'
#        
#        dictdsregion = {}
#        
#        for dsname in dschecklist:
#            
#            if 'match' in dsname:
#                dsnametmp = dsname.split("_")[0]
#                dl = "/".join(dictfiles[dsnametmp]['filelist'][0].split("/")[:-1]+[""])
#                filetag = "_".join(dictfiles[dsnametmp]['filelist'][0].split("/")[-1].split("_")[:-2])+"_match"
#                
#            else:
#                dl = "/".join(dictfiles[dsname]['filelist'][0].split("/")[:-1]+[""])
#                filetag = "_".join(dictfiles[dsname]['filelist'][0].split("/")[-1].split("_")[:-2])
#            
#            if (os.path.exists(dl+"_".join([datasettag,filetag])+'.nc')):
#                
#                # Load from file
#                dictdsregion[dsname] = xr.open_dataset(dl+"_".join([datasettag,filetag])+'.nc')
#            
#            else:
#                dslist = []
#                for regionname in regionlist:
#                # Filter out the points that do not fall within the map we're making
#                
#                    dslist.append(sumplt.dataset_pickpolygon_split(dsdict[dsname],regionname,dir_shape=datadir + 'biomes/'))
#                    
#                ds_tosave = xr.concat(dslist,dim=pd.Index(regionlist, name='regions'))
#                dictdsregion[dsname] = ds_tosave
#                
#                # Save out to netCDF
#                if datasave is True:
#                    
#                    ds_tosave.to_netcdf(path=dl+"_".join([datasettag,filetag])+'.nc')
#            
#        #%% Summary of parameters by region
#        
#        
#        parnamefix = [pn for pn in lparnames if pn != 'IWUE: GPP*VPD/ET: gC/kgH2o *hPa*']
#        
#        # Pick standardization as median of all points of variable
#        normlist = []
#        for pn in parnamefix:
#            #normlist.append(np.nanmedian(dsdict[dschecklist[0]][pn].values.flatten()))
#            normlist.append(np.nanmedian(np.nanmedian(dsdict[dschecklist[0]][pn].values,axis=2).flatten()))
#            
#        d0 = np.array(normlist)  
#        
#        
#        
#        fig, axarr = plt.subplots(int(len(regionlist)/2),2,figsize=(20,20),sharex=True,sharey=True)
#        for ax,rg in zip(axarr.flatten(),regionlist):
#            
#            
#    
#            
#            for idx, fp in enumerate(dschecklist):
#                boxds = dictdsregion[fp].sel(regions=rg)
#                
#                datalist = []
#                for pn in parnamefix:
#                    tmp = boxds[pn].values.flatten()
#                    datalist.append(tmp[~np.isnan(tmp)])
#                    
#                dataarray = np.stack(datalist,axis=1)
#                    
#                
#                
#                
#                d75 = np.percentile(dataarray,75,axis=0)
#                d50 = np.percentile(dataarray,50,axis=0)
#                d25 = np.percentile(dataarray,25,axis=0)
#                
#                #if idx == 0:
#                    #d0 = np.median(dataarray,axis=0)
#            
#                ax.fill_between(np.arange(len(parnamefix)),d25/d0,d75/d0, alpha=0.4)
#                
#                
#                ax.plot(np.arange(len(parnamefix)),d50/d0,'o-',label=fp)
#            ax.set_title(rg,fontsize=10)
#            ax.legend()
#            #ax.set_xlabel('parameters',fontsize=8)
#            #ax.set_ylabel('xy std of median',fontsize=8)
#            ax.tick_params(labelsize = 8.0)
#            
#        plt.xticks(np.arange(len(parnamefix)))
#        
#        
#        axes = plt.gca()
#    
#        axes.set_ylim([0,8])
#        plt.tight_layout()
#        if plotsave is True:    
#            plotname = 'regionpar'
#                            
#            savetag = "_".join([plotname]+dschecklist)
#            plt.savefig(plotdir + savetag)
#            #plt.savefig(plotdir + plotname + savetag +'.eps')
#            plt.close("all")
#    
#        
#        #%% Analysis of Parameter sets by region
#        
#        
#    
#    
#        dictregion = boxregions_abjpl()
#    
#                            
#        regionname = 'SH South America'
#        lt, ln = dictregion[regionname]['lat'], dictregion[regionname]['lon']
#        
#        fig, axarr = plt.subplots(8,4,figsize=(14,28))
#        axarrf = axarr.flatten()
#        
#        for ap,ky in zip([.7,.3],['GCRUN0','GCRUN1']):
#        
#            dsregion = dsdict[ky].copy().sel(lat=slice(lt[0],lt[1]),lon=slice(ln[0],ln[1]))
#            
#            varlist = []
#            for var in dsregion.data_vars:
#                varlist.append(dsregion[var].values)
#                
#                
#            vararray = np.stack(varlist)
#            
#            varregion = vararray.reshape(vararray.shape[0],vararray.shape[1]*vararray.shape[2]*vararray.shape[3])
#            vrnan = np.isnan(varregion.mean(axis=0))
#            varregionnn = varregion[:,~vrnan]
#            
#    
#            for idx, (ax, var) in enumerate(zip(axarrf,dsregion.data_vars)):
#                
#                ax.hist(varregionnn[idx,:],bins=40,alpha=ap,label=ky,normed=True)
#                ax.set_title(regionname + ' : ' + var)
#                
#        axarrf[0].legend()
#        plt.tight_layout()
#        
#        if plotsave is True:    
#                plotname = 'dist_regionensemble_'
#                                
#                savetag = "_".join([regionname,'GCRUN0','GCRUN1'])
#                plt.savefig(plotdir + plotname + savetag)
#                #plt.savefig(plotdir + plotname + savetag +'.eps')
#                plt.close("all")
#        
#        #%% Workingon bring in a distance measurement between sets of parameters into 
#        # the analysis to show distance of median on a map for all parameters
#        # Similar to a correlation
#        
#        if False:
#        
#            from numpy import linalg as LA
#            
#            dist = LA.norm(a-b)
#            
#            m = np.arange(8).reshape(2,2,2)
#            LA.norm(m, axis=(1,2))
#        
#            LA.norm(m[0, :, :]), LA.norm(m[1, :, :])
#    
        #%% Box plot
        # https://stackoverflow.com/questions/43612687/python-matplotlib-box-plot-two-data-sets-side-by-side
