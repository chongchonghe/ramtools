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
import re
import yt
from .yt_field_descrs import FIELDS
from .utilities import my_yt_load, get_sink_info_from_movie1, get_times_from_movie1
from . import utilities as util

from . import units
# from .utils import units, center, utilities, tools

RAM_DIR = ".."

class NoSinkParticle(Exception):
    """exception: no sink particle found in an output or a .csv file."""

    pass

class Ramses():

    def __init__(self, jobdir=None, jobid=None, ram_dir=None, is_use_jobid=False):
        """Initiate an Ramses instance for a given job. The job can be
        specified by either (1) a jobdir, the path of the job directory, or (2)
        a jobid with a ram_dir (default to RAM_DIR), which is equivalent to
        setting jobdir = f"{ram_dir}/Job{jobid}".

        Args:
            jobdir (str): relative/absolute path to the job directory
            jobid (str): jobid, prepend with 'Job' to get the job directory
            ram_dir (str): the base directory for the RAMSES jobs
            is_use_jobid: not being used.
        """

        # if not is_use_jobid:
        if jobid is None:
            # assert not bool(re.compile(r'[~0-9]').search(jobdir[0])), \
            assert not jobdir[0].isdigit(), \
                "jobdir should be the directory of the job, not a jobid"
            self.jobPath = jobdir
        else:
            if ram_dir is None:
                ram_dir = RAM_DIR
            self.jobPath = f"{ram_dir}/Job{jobid}"
        self.get_units()
        self.ds1 = None
        self.tRelax = self.get_trelax()

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
        return

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
        data = np.loadtxt(fp, delimiter=',', usecols = (0,1,2,3,4,5,6,7,8))
        if data.ndim == 1:
            data = np.array([data])
        return data[:, 1:8]
    def get_sink_acc_rate(self, outputID):
        fp = self.get_sink_path(outputID)
        if not os.path.isfile(fp):
            raise FileNotFoundError
        with open(fp, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                raise NoSinkParticle
        data = np.loadtxt(fp, delimiter=',', usecols = (-2))
        if data.ndim == 1:
            data = np.array([data])
        return data[:, :]
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

    def get_sink_masses(self, outputID=None):
        """ Get sink masses in M_sun of one output (not shifted)
        Replacing get_sink_mass

        Return
        ------
            sink mass: 1-D Array
        Exceptions
        ----------
            NoSinkParticle
            FileNotFoundError

        """

        # particle = self.get_sink_particles(outputID)
        # sink_mass = particle[:, 0]
        # sink_mass *= self.unit_m / units.Msun
        # return sink_mass
        return self.get_sink_particles(outputID)[:, 0] * self.unit_m / units.Msun

    def get_trelax(self):
        nml = f90nml.read(self.jobPath + "/run.sink.nml")
        try:
            t_relax = nml['POISSON_PARAMS']['trelax'] * self.unit_t / \
                          units.Myr
        except KeyError:
            t_relax = 0
        return t_relax

    def get_sink_info_from_movie1(self, sinkid):
        """ Return a list of time and a list of sink masses

        Args
            sinkid: interger, starting from 0

        Return
        ------
        Return a dictionary containing the following keys 'out', 't', 'm', 'x',
        'y', 'z', 'vx', 'vy', 'vz'. All in code units.

        """

        return get_sink_info_from_movie1(f"{self.jobPath}/movie1", sinkid)

    def get_times_from_movie1(self):
        return get_times_from_movie1(f"{self.jobPath}/movie1")

    def get_time(self, outputID, readinfo=False):
        """ Get the time in Myr (not substracting t_relax) of data_id.
        """

        fname = "{0}/output_{1:05d}/info_{1:05d}.txt".format(
            self.jobPath, outputID)
        if not os.path.isfile(fname):
            return -1
        # t = yt.load(fname)['time']
        t = util.read_quant_from_ramses_info(fname, 'time')
        return t * self.unit_t / units.Myr

    def get_first_output(self):
        for i in range(1, 100):
            fname = "{0}/output_{1:05d}/info_{1:05d}.txt".format(
                self.jobPath, i)
            if os.path.isfile(fname):
                return i
        return 0

    def get_all_outputs(self):
        outs = []
        for i in range(1, 1000):
            fname = "{0}/output_{1:05d}/info_{1:05d}.txt".format(
                self.jobPath, i)
            if os.path.isfile(fname):
                outs.append(i)
        return outs

    def get_last_output_backward(self):
        for i in range(99, 0, -1):
            fname = "{0}/output_{1:05d}/info_{1:05d}.txt".format(
                self.jobPath, i)
            if os.path.isfile(fname):
                return i
        return 0

    def find_formation_pos_from_movie(self, sinkid, msg='off'):
        """Find the formation location of the sink particles with given sinkid.
        This is done by finding the sink_*.txt file where the sink first form
        and stops accreting.

        Args:
            sinkid: integer (starting from 1)

        Return:
            pos: (n by 3 array) position of sink particle at formation

        """

        # params
        # thresh in fractional growth
        # thresh = 0.00001
        thresh = np.inf         # no thresh

        mpath = "{}/movie1".format(self.jobPath)
        # if not os.path.isfile("{}/sink_00001.txt".format(mpath)):
        #     print("Movies files not found. Quitting")
        #     return
        i = 0
        count = 0
        #par = -1. * np.ones(14) * np.inf
        m0 = -1.0
        for i in range(10000):
            fn = "{}/sink_{:05d}.txt".format(mpath, i)
            if not os.path.isfile(fn):
                continue
            with open(fn, 'r') as f:
                if not os.fstat(f.fileno()).st_size:
                    continue
                num_lines = len(f.readlines())
                #num_lines = 0
                #for line in f:
                #    num_lines += 1
                if num_lines >= sinkid:
                    pars_new = np.loadtxt(fn, delimiter=',')
                    if num_lines >= 2:
                        parnew = pars_new[sinkid-1, :]
                    else:
                        parnew = pars_new
                    mnew = parnew[1] * self.unit_m / units.Msun
                    if mnew - m0 < thresh:
                        # good, particle not growing, print the previous one
                        #assert par is not None
                        #par[1] *= self.unit_m / units.Msun
                        if msg is True or msg == 'on':
                            print("Found a non-growing sink at id {}, after {}"\
                                " steps after formation".format(sinkid, count-1))
                        return parnew[2:5]
                    #par = parnew
                    m0 = mnew

    def overplot_sink(self, p, out, plot_args={}):
        """Over plot sink particles (as green cross) on top of a YT
        slice/project plot"""

        _plot_args = {'color': 'g'}
        for i in plot_args:
            _plot_args[i] = plot_args[i]
        try:
            poss = self.get_sink_positions(out) / self.boxlen
        except NoSinkParticle or FileNotFoundError:
            poss = []
        for pos in poss:
            p.annotate_marker(pos, coord_system='data', plot_args=_plot_args)

    def overplot_sink_with_id(self, plot, out, center, radius, is_id=True,
                              colors=cm.Greens, withedge=False):
        """
        Args:
            plot (yt plot)
            out (int)
            center (tuple or list): the center of the plot in boxlen units
            radius (double): radius in boxlen units
        """

        try:
            sposs = self.get_sink_positions(out) / self.boxlen
            masses = self.get_sink_masses(out)
            indices = np.arange(len(masses))
            _is_inside = np.max(np.abs(sposs - center), axis=1) < radius
        except NoSinkParticle:
            sposs = []
            masses = []
        lim = [1e-2, 1e2]
        #colors = cm.Greens
        for i in range(len(masses)):
            if not _is_inside[i]:
                continue
            m, pos, indice = masses[i], sposs[i], indices[i]
            mass_scaled = (np.log10(m) - np.log10(lim[0]))/np.log10(lim[1]/lim[0])
            plot.annotate_marker(pos, 'o', coord_system='data',
                                 plot_args={'color': colors(mass_scaled),
                                            's':20, 'zorder':i+10,
                                            'linewidths': 0.8,
                                            'edgecolors': 'k' if withedge else 'face'})
            if is_id:
                plot.annotate_text(pos, str(indices[i]), coord_system='data',
                                   text_args={'color': 'k', 'va': 'center', 'ha': 'center',
                                             'size': 8, 'zorder': i+10+0.5})



class RamsesJob(Ramses):

    def __init__(self, jobid):
        return

