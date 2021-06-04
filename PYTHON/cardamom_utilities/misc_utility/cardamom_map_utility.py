#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 21:30:41 2020

@author: gquetin

From cmip5biascorrect_testpoints.py (duplicated there until it is deleted)

Return the rowcol designation for geosCHEM/CARDAMOM file naming

function for switching between 180 and 0-360


"""

import scipy.io as sio
import numpy as np
import os
import sys
import imp

#imac: '/Users/gregoryquetin'
#sherlock: '/home/users/gquetin'
MACHDIR = os.path.expanduser("~")
sys.path.append(MACHDIR + '/repos/scripts/Python/Projects/J5') 

def latlon_to_rowcol(lat,lon):
    
    # Vectors for code in file name
    latspace = np.linspace(-90,90,46)
    lonspace = np.linspace(-180,180,73)
    
    latidx = np.argwhere(latspace == lat)
    lonidx = np.argwhere(lonspace == lon)
    
    if len(latidx) == 0 or len(lonidx) == 0:
        rowcol = []
        print('No point at lat = ' + str(lat) + ' or lon = ' + str(lon))
    
    else:
        latnum = latidx[0][0] + 1
        lonnum = lonidx[0][0] + 1
        
        if latnum < 10:
            latstr = '0'+ str(latnum)
        else:
            latstr = str(latnum)
            
        if lonnum < 10:
            lonstr = '0' + str(lonnum)
        else:
            lonstr = str(lonnum)
        
        rowcol = latstr+lonstr
        
    return rowcol




def rowcol_to_latlon(checkptlist,plotcheck=False):
    """
    Convert GEOCHEM rows and columns to lats and longs
    """
    
    # Vectors for code in file name
    latspace = np.linspace(-90,90,46)
    lonspace = np.linspace(-180,180,73)
    
    latpts = []
    lonpts = []
    for ckpt in checkptlist:
        
        latpts.append(latspace[int(ckpt[0:2])-1])
        lonpts.append(lonspace[int(ckpt[2:])-1])
    
    latlonzip = list(zip(latpts,lonpts))
    
#    if plotcheck==True:
#        import matplotlib
#        matplotlib.use("Agg") 
#        import matplotlib.pyplot as plt
#        import mpl_toolkits.basemap as bm #This has a problem with 'pyproj'
#        from mpl_toolkits.basemap import Basemap
#        #create instance of Base map centered at 180
#        centerlon = 0;
#        m = Basemap(projection='moll',lon_0=centerlon,lat_0 = 0,resolution='c') #'l' is for low res.,'c' is for corse
#        
#        
#        #create figure and axes handles
#        fig = plt.figure()
#        ax = fig.add_axes([0.05,0.05,0.9,0.9])
#        
#        #plot data
#        for k,(lat,lon) in enumerate(latlonzip):
#            #m.scatter(lon,lat,color = 'r',s=100,latlon=True)
#            x,y = m(lon,lat)
#            plt.annotate(checkptlist[k], xy=(x, y), xytext=(x,y),
#                         color='r',fontsize=10)
#        
#        #add lon and lat lines
#        m.drawparallels(np.arange(-90.,90.,30.),labels=[1,1,0,0]) #labels = [left,right,top,bottom]
#        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0])
#        
#        #draw map boundary
#        m.drawcoastlines(linewidth=1.25)
#        m.drawmapboundary(fill_color='1') #make map background white.To make it black, 
#        
#        #draw country borders
#        m.drawcountries()
#        
#        #add title
#        ax.set_title('lat lon pts')
#        
#        #display plot on screen
#        plt.show()
    
    
    return latlonzip


def lon180_to_lon360(lat,lon):

    dslat = lat
    
    # 
    if lon < 0:
        dslon = lon + 360
    else:
        dslon = lon
    
    return dslat,dslon