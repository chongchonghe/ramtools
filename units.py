#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Define some astronomical constants in cgs units

@author: chongchonghe
"""

import yt.units as U

# constants in cgs units
Msun = float(U.Msun.in_cgs())
pc = float(U.pc.in_cgs())
AU = float(U.AU.in_cgs())
Myr = float(U.Myr.in_cgs())
G = float(U.G.in_cgs())
mH = float(U.mass_hydrogen_cgs)
