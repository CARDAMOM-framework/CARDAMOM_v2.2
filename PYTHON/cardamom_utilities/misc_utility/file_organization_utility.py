#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 22:03:50 2020

@author: gquetin

This script is to contain specific file organization for running CARDAMOM

"""

import os

def folderdeff(homedir):
    
    if homedir == '/Users/gregoryquetin':
        datadir = homedir+'/Google Drive/DATA/'
        dataoutdir = homedir+'/Google Drive/DATA/Analysis/'
        #testdir = homedir+'/Google Drive/notadirectory'
        
    elif homedir == '/home/users/gquetin':
        datadir ='/scratch/users/gquetin/DATA/'
        dataoutdir = '/scratch/users/gquetin/DATA/Analysis/'
        
    elif homedir == '/Users/gquetin':
        datadir = homedir+'/Google Drive/DATA/'
        dataoutdir = homedir+'/Google Drive/DATA/Analysis/'
        
    dirlist = [datadir,dataoutdir]#,testdir]
    
    for dl in dirlist:
        if not os.path.exists(dl):
            os.makedirs(dl)
            
    return datadir, dataoutdir