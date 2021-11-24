# -*- coding: utf-8 -*-
"""
Defines the Ramses class, the core of ramtools.

Attributes:
    RAM_DIR (str): Global variable defining the path to RAMSES jobs
"""

import os
import warnings
import json
import numpy as np
from matplotlib.pyplot import cm
from glob import glob
import f90nml
import yt
from . import utilities as util
from . import units

RAM_DIR = None

class Ramses():
    """This is the core of ramtools. Most of the times, start using ramtools  
by claiming this class. 

Examples
--------
>>> import ramtools as rt
>>> ram = rt.Ramses("../tests/Job_test")
>>> par = ram.get_sink_particles(1)
>>> mass = par[:, 0]
>>> print("Total mass =", mass.sum())
    
>>> rt.set_RAM_DIR("../tests")
>>> ram = rt.Ramses(jobid="1")

"""

    def __init__(self, jobdir=None, jobid=None, ram_dir=None):
        """Declare a Ramses instance for a given job. The job can be
        specified by one of the following options:

        1. r = Ramses(`job_dir`), where `job_dir` is the path to the job directory. 
        2. r = Ramses(jobid = `_id`, ram_dir = `ram_dir`), where `ram_dir` is the directory where all the RAMSES jobs are stored, and `_id` is the string after "Job". This is equivalent to Ramses(`ram_dir/Job_id`)
        3. ramtools.set_RAM_DIR(`ram_dir`); r = Ramses(jobid = `_id`).

        Args:
            jobdir : str
                relative/absolute path to the job directory
            jobid : str
                jobid, postfix to 'Job' to get the job directory name (in `ram_dir`)
            ram_dir : str
                the base directory for the RAMSES jobs
        """

        if jobdir is not None:
            self.jobPath = jobdir
        elif jobid is None:
            warnings("Need to specify either jobdir or jobid")
            return
        else:
            if ram_dir is None:
                if RAM_DIR is None:
                    warnings("Need to define ram_dir. Either specify"
                    "ram_dir as a parameter, or define a global RAM_DIR"
                    "by calling ramtools.set_RAM_DIR(`ram_dir`)")
                    return
                else:
                    ram_dir = RAM_DIR
            self.jobPath = f"{ram_dir}/Job{jobid}"
        self.ds1 = None
        if self.get_units() == 0:
            self.tRelax = self.get_trelax()
        self.ramses_dir = os.path.dirname(self.jobPath if self.jobPath[-1] !=
                                          '/' else self.jobPath[:-1])
        self.data_dir = os.path.join(self.ramses_dir, "h5_data")

    def get_info_path(self, out):
        """Return the path to info_out.txt"""
        return "{0}/output_{1:05d}/info_{1:05d}.txt".format(
            self.jobPath, out)

    def load_ds(self, out):
        """Return a yt.load instance of the frame `out`"""
        return yt.load(self.get_info_path(out), fields=util.FIELDS)

    def get_ds(self):
        """Load the first output out there. This is necessary to get the units
        and other things"""
        for i in range(1, 100):
            if os.path.isfile(self.get_info_path(i)):
                self.ds1 = self.load_ds(i)
                return 0        # success
        return 1                # fail to find an output

    def get_units(self):
        """Define the following units for this job: 

        1. unit_l: this is the boxlen (in cm), which equals to unit_l_code * boxlen
        2. unit_l_code: this is the actually length unit (in cm)
        3. unit_d: dnesity unit in cgs
        4. unit_t: time unit in cgs
        5. unit_v: velocity unit in cgs
        6. unit_m: mass unit in cgs
        7. unit_m_in_Msun: mass unit in Msun
        8. kin_ene_in_cgs: kinectic energy in cgs
        """
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
        self.unit_B = util.get_unit_B(self.ds1)
        # self.n_colden_H = np.double(self.unitDen2Hcc * self.unit_l)
        return 0

    def get_sink_path(self, out):
        """Return the path to sink_*.csv """
        sinkFile = "{0}/output_{1:05d}/sink_{1:05d}.csv".format(
            self.jobPath, out)
        return sinkFile

    def get_sink_particles(self, out):
        """
        Args:
            out (int): the output frame

        Return: 
            particles (array): (n, 7) array containing the sink
                parameters, where n is the number of particles and the 7
                columns are: [m, x, y, z, vx, vy, vx], all in code units

        Raise:
            FileNotFoundError
            NoSinkParticle

        """

        fp = self.get_sink_path(out)
        if not os.path.isfile(fp):
            raise FileNotFoundError
        with open(fp, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                raise NoSinkParticle
        data = np.loadtxt(fp, delimiter=',', usecols = (0,1,2,3,4,5,6,7,8))
        if data.ndim == 1:
            data = np.array([data])
        return data[:, 1:8]

    def get_sink_positions(self, out):
        """
        
        Args:
            out (int): the output frame

        Return:
            array: an array with shape (n, 3) containing the particle positions 
                in code unit

        Raise:
            FileNotFoundError
            NoSinkParticle

        """

        return self.get_sink_particles(out)[:, 1:4]

    def get_sink_masses(self, outputID=None):
        """ Get sink masses in M_sun of one output (not shifted)
        Replacing get_sink_mass

        Args:
            out (int): the output frame

        Return:
            array: an array of shape (n, ) containing the particle masses in
                **solar mass**.

        Raise:
            FileNotFoundError
            NoSinkParticle


        """

        return self.get_sink_particles(outputID)[:, 0] * self.unit_m / units.Msun

    def get_sink_acc_rate(self, out):
        """Get the sink accretion rate. 
        
        .. warning:: DO NOT USE THIS. The results are not trustable. This column of data from RAMSES outputs seems to be problematic

        """

        fp = self.get_sink_path(out)
        if not os.path.isfile(fp):
            raise FileNotFoundError
        with open(fp, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                raise NoSinkParticle
        data = np.loadtxt(fp, delimiter=',', usecols = (-2))
        if data.ndim == 1:
            data = np.array([data])
        return data[:, :]

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

        Args:
            sinkid (int): sink id, starting from 0

        Returns:
            int: a dictionary containing the following keys: 'out', 't', 'm', 'x',
                'y', 'z', 'vx', 'vy', 'vz'. All in code units.

        """

        return util.get_sink_info_from_movie1(f"{self.jobPath}/movie1", sinkid)

    def get_times_from_movie1(self):
        return util.get_times_from_movie1(f"{self.jobPath}/movie1")

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
            sinkid (int): (starting from 1)
            msg (bool or string): print message if msg is True or 'on'

        Return:
            array: n by 3 array. Position of sink particle at formation

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
        """Over plot sink particles (as green crosses) on top of a YT
        slice/project plot
        
        Args:
            p (yt.sliceplot or yt.projectplot): the plot to overplot on
            out (int): the output frame
            plot_args (dict): 

        """

        _plot_args = {'color': 'g'}
        plot_args.update(_plot_args)
        try:
            poss = self.get_sink_positions(out) / self.boxlen
        except NoSinkParticle or FileNotFoundError:
            poss = []
        for pos in poss:
            p.annotate_marker(pos, coord_system='data', plot_args=plot_args)

    def overplot_sink_with_id(self, plot, out, center, radius, is_id=True,
                              colors=cm.Greens, withedge=False):
        """
        Args:
            plot (yt plot): the plot to overplot on
            out (int): the output frame
            center (tuple or list): the center of the plot in boxlen units
            radius (float): the radius (half width) of the box around
                the center defining the domain of interest, in boxlen units

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


class NoSinkParticle(Exception):
    """exception: no sink particle found in an output or a .csv file."""
    pass


def set_RAM_DIR(ram_dir):
    """Set the global variable RAM_DIR. Check Ramses.__init__ for its usage."""
    
    global RAM_DIR
    RAM_DIR = ram_dir