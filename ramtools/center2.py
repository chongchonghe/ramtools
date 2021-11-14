#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Use center.py. Jan 21, 2019
#Replacing center.py

@author: chongchonghe
"""

# constants
# HYDRO_FRAC = 0.76
sigma0 = 6.304e-18  # TODO: add sigma0 for HeI and HeII
mH = 1.6733e-24  # g
mHe = 6.6465e-24
# nH2nHe = ((1 - HYDRO_FRAC) / mHe) / (HYDRO_FRAC / mH)
eHI = 13.6
eHeI = 24.6
eHeII = 54.4
atoms = ['HI', 'HeI', 'HeII']
unit_lum = 1e44 # 1e44 ergs/s
qscale = 1.0 / unit_lum

LOG_TAO_LIM = [-1.0, 2.0] # TODO fix logtauLim

jobnames = {
    '4.01':      'XSF',
    '4.0.1':     'SF', 
    '4.2.1':     'MF', 
    '4.3.1':     'LF', 
    # '4.3':       'LF',
    '4.4':       'XLF',
    # '4.4.1':     'XLF',
    '2.01':      'XSC',
    '2.0':       'SC', 
    '2.2.2':     'MC', 
    '2.3':       'LC', 
    '3.001':     'XXSV', 
    '3.01':      'XSV', 
    '3.0':       'SV', 
    '3.2.2':     'MV', 
    '3.3.2':     'LV', 
    # '2.3.lm':    'LC',
    # '2.3.lm':    'None', 
    # '2.3.lm2':   'None', 
    }

grid_fig_size = [2.4 * 4.4, 1.8 * 4.4 * 3/4]

plot_range_in_tff = [0, 7.0] # [0, plotrangetff * tff]


# v2: Corrected Kroupa normalization; set mass_shift to 0.4, 2018/12/13
global_version = '_v2'

jobids = {'4':['4.01', '4.0.1', '4.2.1', '4.3.1', '4.4'],
          '2':['2.01', '2.0',   '2.2.2', '2.3'],
          # '3':['3.01', '3.0',   '3.2.2', '3.3.2']}
          '3':['3.001','3.01',  '3.0',   '3.2.2', '3.3.2']}
