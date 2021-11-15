#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Attributes:

    FIELDS (list): a list of RAMSES fields. This is different from the default 
        fields in the public version of RAMSES

"""


FIELDS = ["Density",
        "x-velocity", "y-velocity", "z-velocity",
        "x-Bfield-left", "y-Bfield-left", "z-Bfield-left",
        "x-Bfield-right", "y-Bfield-right", "z-Bfield-right",
        "Pressure",
        "xHII", "xHeII", "xHeIII"]

# from yt.frontends.ramses.data_structures import *
# RAMSESIndex.fluid_field_list = FIELDS