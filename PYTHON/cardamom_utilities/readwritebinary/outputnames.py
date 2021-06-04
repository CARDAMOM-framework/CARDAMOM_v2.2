#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 22:00:34 2020

@author: gquetin

Output names and lists for CARDAMOM

"""

def poolnames(output='longnames'):
    
    # DEFINITIONS
    longnames = ['Labile C',
                 'Foliar C',
                 'Root C',
                 'Wood C',
                 'Fine litter C',
                 'Soil organic matter (SOM) C',
                 'Plant-available H2O'
                 ]
    
    shortnames = ['c_labile',
                  'c_foliar',
                  'c_root',
                  'c_wood',
                  'c_finelitter',
                  'c_som',
                  'h2o_forplant'
                  ]
    
    units = ['gC/m2',
             'gC/m2',
             'gC/m2',
             'gC/m2',
             'gC/m2',
             'gC/m2',
             'kgH20/m2'
             ]
    
    if output == 'longnames':
        return longnames
    elif output == 'shortnames':
        return shortnames
    elif output == 'dict':
        dictname = {}
        for sn, ln in zip(shortnames,longnames):
            dictname[sn]=ln
        
        return dictname
    
    elif output == 'dictunits':
        dictunit = {}
        for sn, unit in zip(shortnames,units):
            dictunit[sn]=unit
        
        return dictunit

def fluxnames(output='longnames'):
    
    longnames = ['GPP',
              'temprate (decomposition factor)*',
              'respiration_auto (GPP -> autotrophic respiration)',
              'leaf_production (GPP -> leaves)',
              'labile_production (GPP -> labile)',
              'root_production (GPP -> roots)',
              'wood_production (GPP -> wood)',
              'labile_release (Labile -> foliar)',
              'leaffall_factor (leaf senescence factor)*',
              'leaflitter_production (leaf -> litter)',
              'woodlitter_production  (wood -> soil organic C)',
              'rootlitter_production  (root -> litter)',     
              'respiration_het_litter (litter respiration)',
              'respiration_het_som (SOM respiration)',
              'litter2som (litter -> SOM)',
              'labrelease_factor(leaf onset factor)*',
              'Fires (total fire emissions)',
              'Fires (fire C fluxes from Labile)',
              'Fires (fire C fluxes from Foliar)',
              'Fires (fire C fluxes from Root)',
              'Fires (fire C fluxes from Wood)',
              'Fires (fire C fluxes from Fine litter)',
              'Fires (fire C fluxes from Soil organic matter)',
              'Fires (C pool transfers, to litter and SOM pools)',
              'Fires (C pool transfers, to litter and SOM pools)',
              'Fires (C pool transfers, to litter and SOM pools)',
              'Fires (C pool transfers, to litter and SOM pools)',
              'Fires (C pool transfers, to litter and SOM pools)',
              'ET',
              'Runoff']
    
    shortnames = ['gppflux',
                  'decf_tempr',
                  'gpp_to_autoresp',
                  'gpp_to_leaf',
                  'gpp_to_labile',
                  'gpp_to_root',
                  'gpp_to_wood',
                  'labile_to_foliar',
                  'leaf_fall',
                  'leaf_to_litter',
                  'wood_to_soilc',
                  'root_to_litter',
                  'hetresp_litter',
                  'hetresp_som',
                  'litter_to_som',
                  'leaf_onset',
                  'fire_em_total',
                  'fire_em_labile',
                  'fire_em_foliar',
                  'fire_em_root',
                  'fire_em_wood',
                  'fire_em_litter',
                  'fire_em_som',
                  'fire1_to_littersom',
                  'fire2_to_littersom',
                  'fire3_to_littersom',
                  'fire4_to_littersom',
                  'fire5_to_littersom',
                  'et',
                  'runoff'
                  ]
    
    units = ['gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'gC/m2/day',
             'kgH20/m2/day',
             'kgH20/m2/day'
             ]
    
    if output == 'longnames':
        return longnames
        
    elif output == 'shortnames':
        return shortnames

    elif output == 'dict':
    
        dictname = {}
        for sn, ln in zip(shortnames,longnames):
            dictname[sn]=ln
        
        return dictname
    
    elif output == 'dictunits':
        dictunit = {}
        for sn, unit in zip(shortnames,units):
            dictunit[sn]=unit
        
        return dictunit

def derivedfluxnames(output='shortnames'):
    
    longnames = ['Leaf Area Index',
                 'Gross Primary Productivity',
                 'Net Biosphere Exchange',
                 'Net Ecosystem Exchange'
                 ]
    
    
    
    shortnames = ['lai',
                  'gpp',
                  'nbe',
                  'nee']
    
    units = 'gC/m2/day'
    
    if output == 'longnames':
        return longnames
        
    elif output == 'shortnames':
        return shortnames

    elif output == 'dict':
    
        dictname = {}
        for sn, ln in zip(shortnames,longnames):
            dictname[sn]=ln
        
        return dictname
    
    elif output == 'units':
        return units
    
# %%

if __name__ == '__main__':
    
    print('No main file functions currently, check cardamom_output2netcdf.py')