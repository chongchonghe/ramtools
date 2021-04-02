#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""center.py
Define gloable variables for the RAMSES project.

Author: chongchonghe
"""

from __future__ import division
import os

MASS_SHIFT = 0.4
den_c_to_m = 0.194556           # central density to mean density
PLOT_RANGE_IN_TFF = [0, 7.0] # [0, plotrangetff * tff]

JOBNAMES = {
    '4.01': 'XS-F', '4.0.1': 'S-F', '4.2.1': 'M-F', '4.3.1': 'L-F', '4.4': 'XL-F',
    '2.01': 'XS-C', '2.0': 'S-C', '2.2.2': 'M-C', '2.3': 'L-C',
    '2.3.lm2': 'L-C-lm', '2.3.lm': 'L-C-xlm',
    '3.001': 'XXS-VC', '3.01': 'XS-VC', '3.0': 'S-VC', '3.2.2': 'M-VC', '3.3.2': 'L-VC',
    '4.0': 'S-F', '4.2': 'M-F', '4.3': 'L-F', '4.4.1': 'XL-F',
    '3.2.2.t': 'M-VC',
    '1.2.5': 'try',
}

