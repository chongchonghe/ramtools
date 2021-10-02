#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Define the following constants in cgs units:
    - all constants in con_list
    - all constants in units_list
    - m_H

@author: chongchonghe
"""

# import yt.units as U

# # constants in cgs units
# Msun = float(U.Msun.in_cgs())
# pc = float(U.pc.in_cgs())
# AU = float(U.AU.in_cgs())
# Myr = float(U.Myr.in_cgs())
# G = float(U.G.in_cgs())
# mH = float(U.mass_hydrogen_cgs)
# yr = float(U.yr.in_cgs())

from astropy import units as U
from astropy import constants as C

con_list = ['G', 'N_A', 'R', 'Ryd', 'a0', 'alpha', 'atm', 'b_wien', 'c', 'g0',
           'h', 'hbar', 'k_B', 'm_e', 'm_n', 'm_p',
           'sigma_T', 'sigma_sb', 'u', 'GM_earth', 'GM_jup', 'GM_sun',
           'L_bol0', 'L_sun', 'M_earth', 'M_jup', 'M_sun', 'R_earth', 'R_jup',
           'R_sun', 'au', 'kpc', 'pc']
# Some EM constants is unable to load because different cgs system
# has different values

con_list.sort(key=lambda y: y.lower())
for con in con_list:
    locals()[con] = eval("C.{}.cgs.value".format(con))

units_list = ['yr', 'Myr', 'Gyr']
for _var in units_list:
    locals()[_var] = eval("U.{}.cgs._scale".format(_var))

m_H = 1.6733e-24
mH = m_H
Msun = M_sun
AU = au
