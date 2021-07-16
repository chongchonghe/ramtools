#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate luminosity histograms of sink particles in a simulation, for each outputs.

Example
-------
from ramtools.luminosity import lumhist
lumhist(jobdir="../../2017-RAMSES/Job2.2.2/", recipe="frig_he", plotdir="./trash", xlim=[38, 52])

"""

#-----------------------------------------------------------------------------
#    Author: Chong-Chong He
#    Data: 2021-07-02
#-----------------------------------------------------------------------------

from __future__ import division, print_function
import os
import numpy as np
import matplotlib.pyplot as plt

from . import ramses
from . import radiation_module as rm
from . import units

def lumhist(jobdir, recipe, plotdir, xlim=[38, None], ylim=None):
    """
    Args:
        jobidr (str): path to job directory
        recipe (str): one of "frig_he", "sp"
        plotdir (str): directory to save plots
        xlim: xlim of the histogram
        ylim: ylim
    """

    r = ramses.Ramses(jobdir=jobdir)
    os.makedirs(plotdir, exist_ok=True)
    for out in r.get_all_outputs():
        try:
            ms = r.get_sink_masses(out)
        except ramses.NoSinkParticle:
            continue
        lum = rm.luminosity(ms, recipe, r.get_info_path(out))
        # bins = np.logspace(40, 50, 11)
        bins = 'auto'
        plt.hist(np.log10(lum), bins=bins)
        plt.xlabel("log L")
        plt.ylabel("N")
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.savefig(plotdir + f"/output_{out:05d}.pdf")
        plt.close()
