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

def lumhist(jobdir, recipe, plotdir, xlim=[38, None], ylim=None, bins = 'auto'):
    """
    Args:
        jobidr (list of str): paths to job directories
        recipe (str): one of "frig_he", "sp"
        plotdir (str): directory to save plots
        xlim: xlim of the histogram
        ylim: ylim
    """
    print(jobdir, type(jobdir))
    if isinstance(jobdir, str):
        jobdir = [jobdir] #Turns a single directory string into a single item list
    r = ramses.Ramses(jobdir=jobdir[0])
    os.makedirs(plotdir, exist_ok=True)
    colors = ['tab:blue', 'tab:green','tab:red','tab:purple']

    # make L(m) curve
    ms = np.logspace(-1, 3, 100)
    plt.plot(np.log10(ms), np.log10(rm.QVaccaRaw(ms)), 'k--', label="Vacca(M)")
    plt.plot(np.log10(ms), np.log10(rm.QVacca_sp(0.4*ms)), 'k-.',
             label="Vacca_sp(0.4 M)")
    for mtot in np.logspace(0, 4, 5):
        lum = rm.luminosity(ms, recipe, r.get_info_path(1),
                           masstot=mtot)
        plt.plot(np.log10(ms), np.log10(lum),
                 label=r"$M_{{\rm tot}} = {:.0f} M_\odot$".format(mtot))
        plt.gca().set(xlabel="log M [Msun]", ylabel="log L",
                     xlim=[-1, 3])
    plt.legend()
    plt.savefig(plotdir + "/luminosity curve.pdf")
    plt.close()
#    return

    for out in r.get_all_outputs():
        for i, dir in enumerate(jobdir):
            r = ramses.Ramses(jobdir=dir)
            try:
                ms = r.get_sink_masses(out)
            except ramses.NoSinkParticle:
                continue
            lum = rm.luminosity(ms, recipe, r.get_info_path(out))
            # bins = np.logspace(40, 50, 11)
    #        bins = 'auto'
            print(bins)
            plt.hist(np.log10(lum), bins=bins, alpha=0.5, label = f'{dir}  L_tot: {np.sum(lum):.2e}', color = colors[i])
            plt.xlabel("log L")
            plt.ylabel("N")
            plt.title(f"Luminosity histogram at output {out}")
            plt.legend()
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.savefig(plotdir + f"/output_{out:05d}.pdf")
        plt.close()
