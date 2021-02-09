#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools to post-process RAMSES data

Examples:
>>> ram = Ramse("../jobs/Job2")
>>> par = ram.get_sink_particles(1)
>>> mass = par[:, 0]
>>> print("Total mass =", mass.sum())


@author: ChongChong He
"""

import os
import json
import numpy as np
from matplotlib.pyplot import cm
from glob import glob
import f90nml
import yt
from .yt_field_descrs import FIELDS
from .utilities import my_yt_load
from . import units
# from .utils import units, center, utilities, tools

class NoSinkParticle(Exception):
    """exception: no sink particle found in an output or a .csv file."""

    pass


class Ramses():

    def __init__(self, jobdir):
        """
        Args:
            jobdir (str): relative/absolute path to job directory
        """

        self.jobPath = jobdir
        self.get_units()
        self.ds1 = None

    def get_info_path(self, out):
        return "{0}/output_{1:05d}/info_{1:05d}.txt".format(
            self.jobPath, out)

    def load_ds(self, out):
        return yt.load(self.get_info_path(out), fields=FIELDS)

    def get_ds(self):
        for i in range(1, 100):
            if not os.path.isfile(self.get_info_path(i)):
                continue
            self.ds1 = self.load_ds(i)
            return 0            # success
        return 1                # fail to find an output

    def get_units(self):
        flag = self.get_ds()           # get self.ds1
        if flag:
            print("Failed to load an output.")
            return 1
        self.unit_l = np.double(self.ds1.length_unit.value) # box_len = 1
        self.boxlen = self.ds1['boxlen']
        self.unit_l_code = self.unit_l / self.boxlen
        self.unit_d = np.double(self.ds1.density_unit.value)
        self.unit_t = np.double(self.ds1.time_unit.value)
        self.unit_v = self.unit_l_code / self.unit_t
        self.unit_m = np.double(self.ds1.mass_unit.value)
        self.unit_m_in_Msun = np.double(self.ds1.mass_unit.to('Msun').value)
        # self.unitDen2Hcc = self.ds1.density_unit * center.mfrac_H / units.mH
        self.kin_ene_in_cgs = self.unit_m * self.unit_l_code**2 / self.unit_t**2
        self.pot_ene_in_cgs = units.G * self.unit_m**2 / self.unit_l_code
        # self.n_colden_H = np.double(self.unitDen2Hcc * self.unit_l)

    def get_sink_path(self, outputID):
        sinkFile = "{0}/output_{1:05d}/sink_{1:05d}.csv".format(
            self.jobPath, outputID)
        return sinkFile

    def get_sink_particles(self, outputID):
        """Return a (n, 7) array of sink parameters, where n is the number of
        particles and the 7 columns are: [m, x, y, z, vx, vy, vx], all in code
        units

        Args:
            outputID (int):

        Return:
            particles (n by 7 array): particles parameters

        Raise:
            FileNotFoundError
            NoSinkParticle

        """

        fp = self.get_sink_path(outputID)
        if not os.path.isfile(fp):
            raise FileNotFoundError
        with open(fp, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                raise NoSinkParticle
        data = np.loadtxt(fp, delimiter=',')
        if data.ndim == 1:
            data = np.array([data])
        return data[:, 1:8]

    def get_sink_positions(self, outputID):
        """Read n by 3 array of sink particle positions, where n is the number of
        particles: [x, y, z], all in code units (~ pc)

        Return
        ------
        - An array with ndim=2

        Raise
        -----
        - FileNotFoundError
        - NoSinkParticle

        """

        return self.get_sink_particles(outputID)[:, 1:4]

    def get_trelax(self):
        nml = f90nml.read(self.jobPath + "/run.sink.nml")
        try:
            t_relax = nml['POISSON_PARAMS']['trelax'] * self.unit_t / \
                          units.Myr
        except KeyError:
            t_relax = 0
        return t_relax
