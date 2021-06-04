#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:34:43 2019

@author: gquetin

Store settings for the various .cbr related CARDAMOM outputs

"""

def cbr_par_dictionary():
    
    '''
    Contains a dictionary to convert long names from 'C file' into short names
    good for netcdf.
    '''
    pardict = {'Decomposition rate':('decomp_rate','gC/m2/day'),
             'Fraction of GPP respired':('fgpp_r', 'fraction'),
             'Fraction of (1-fgpp) to foliage':('fgpp_to_foliage','fraction'),
             'Fraction of (1-fgpp) to roots':('fgpp_to_roots','fraction'),
             'Leaf Lifespan':('leaf_lifespan','days'), 
             'TOR wood* - 1% loss per year value':('tor_wood','rate'),
             'TOR roots':('tor_root','rate'),
             'TOR litter':('tor_litter','rate'),
             'TOR SOM':('tor_som','rate'),
             'Temp factor* = Q10 = 1.2-1.6':('q10','rate'),
             'Canopy Efficiency':('canopy_eff','fraction'),
             'Bday':('bday','doy'),
             'Fraction to Clab':('f_to_clab','fraction'),
             'Clab Release period':('clab_rperiod','duration'),
             'Fday':('fday','doy'),
             'Leaf fall period':('leaf_fperiod','duration'),
             'LMCA':('lmca','unknown'),
             'C labile':('c_labile','gC/m2'),
             'C foliar':('c_foliar','gC/m2'),
             'C roots':('c_root','gC/m2'),
             'C_wood':('c_wood','gC/m2'),
             'C litter':('c_litter','gC/m2'),
             'C_som':('c_som','gC/m2'),
             'IWUE: GPP*VPD/ET: gC/kgH2o *hPa':('iwue','hPa gC/kgH2O'),
             'Runoff focal point (~maximum soil storage capacity x 4)':('runoff_fpt','unknown'),
             'Runoff focal point':('runoff_fpt','unknown'),
             '"Wilting point"':('wilt_pt','unknown'),
             '"Bucket at t0"':('bucket_t0','kgH2O/m2'),
             'Foliar biomass CF':('folia_biomass_cf','unknown'),
             '"Ligneous" biomass CF".':('ligneous_biomass_cf','unknown'),
             'DOM CF".':('dom_cf','unknown'),
             'Resilience factor (since transfer to litter is represented as (1-pars[30])) ".':('res_fac','unknown'),
             'Lab pool lifespan':('lab_lifespan','duration'),
             'Moisture factor':('moisture_factor','fraction'),
             'Fire-induced mortality':('fire_mortality','unknown'),
             ' ball-berry slope m [unitless]':('m_stomatal','unknown'),
             ' ball-berry intercept b [mol m^-2 s^-1]':('b_stomatal','unknown'),
             ' leaf boundary layer conductance gb [mol m^-2 s^-1]':('gb','mol m^-2 s^-1'),
             ' Jmax [umol m^-2 s^-1]':('Jmax','umol m^-2 s^-1'),
             ' Vmax [umol m^-2 s^-1]':('Vmax','m^-2 s^-1')#803
             }
    
    
    return pardict


def grouped_parnames(group = None):
    '''
    Manually created groupings of the parameters for analysis. Returns a list 
    of parameters for:
        carbon pools
        water pools
        turn over rates
        life span
        critical parameters
        gpp fractions
        water related
        unsorted
        canopy
    '''
     
    carbon_pools = ['c_labile',
                  'c_foliar',
                  'c_root',
                  'c_wood',
                  'c_litter',
                  'c_som',
                  ]

    water_pools = ['bucket_t0']    
        
    turn_over_rates = ['tor_wood',
                       'tor_root',
                       'tor_litter',
                       'tor_som']
    
    life_span = ['leaf_lifespan',
                 'lab_lifespan']
    
    critical_params = ['decomp_rate',
                       'q10',
                       'canopy_eff']
    
    gpp_fractions = ['fgpp_r',
                     'fgpp_to_foliage',
                     'fgpp_to_roots',
                     'f_to_clab']
    
    water_related = ['wilt_pt',
                     'moisture_factor',
                     'iwue',
                     'runoff_fpt'
                     ]
    
    unsorted = ['bday',
                'clab_rperiod',
                'fday',
                'leaf_fperiod',
                'lmca', # leaf mass carbon per area?
                'ligneous_biomass_cf',
                'dom_cf',
                'fire_mortality']
    
    canopy = ['leaf_lifespan',
              'fday',
              'leaf_fperiod',
              'c_foliar',
              'fgpp_to_foliage',
              'canopy_eff']
    
    groups_dict = {'carbon_pools':carbon_pools,
                   'water_pools':water_pools,
                   'pools': carbon_pools + water_pools,
                   'turn_over_rates':turn_over_rates,
                   'life_span':life_span,
                   'critical_params':critical_params,
                   'gpp_fractions':gpp_fractions,
                   'water_related':water_related,
                   'unsorted':unsorted,
                   'canopy':canopy,
                   'soil':['tbd']} 
    
    if group == None:
        
        output =  list(groups_dict.keys())
        
    else:
        output = groups_dict[group]
        
    return output

