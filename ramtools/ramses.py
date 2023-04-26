# -*- coding: utf-8 -*-
"""
Defines the Ramses class, the core of ramtools.

Attributes:
    RAM_DIR (str): Global variable defining the path to RAMSES jobs
"""

from . import center
from . import plotutils as pltu
from .utilities import my_yt_load, read_zoom_radius, read_zoom_center
from . import utilities
from .ramsesbase import *
import warnings
import f90nml
try:
    from .yt_field_descrs import FIELDS
except:
    FIELDS = None

RAM_DIR = None

class Ramses(RamsesBase):
    """Inherit from RamsesBase for the use of the author's personal use. The
    methods defined in this class may not work for another person's RAMSES
    simulations."""

    def __init__(self, jobdir=None, jobid=None, ram_dir=None, fields=FIELDS):
        """Declare a Ramses instance for a given job. The job can be
        specified by one of the following options:

        1. r = Ramses(`job_dir`), where `job_dir` is the path to the job directory.
        2. r = Ramses(jobid = `_id`, ram_dir = `ram_dir`), where `ram_dir` is the directory where all the RAMSES jobs are stored, and `_id` is the string after "Job". This is equivalent to Ramses(`ram_dir/Job_id`)
        3. ramtools.set_RAM_DIR(`ram_dir`); r = Ramses(jobid = `_id`).

        Args:
            jobdir (str): relative/absolute path to the job directory
            jobid (str): jobid, postfix to 'Job' to get the job directory
                name (in `ram_dir`)
            ram_dir (str): the base directory for the RAMSES jobs

        """
        if jobdir is not None:
            jobPath = jobdir
            jobid = os.path.basename(jobdir)[3:]
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
            jobPath = f"{ram_dir}/Job{jobid}"
        RamsesBase.__init__(self, jobPath, fields=fields)
        self.jobid = jobid
        self.ramses_dir = os.path.dirname(self.jobPath if self.jobPath[-1] !=
                                          '/' else self.jobPath[:-1])
        self.data_dir = os.path.join(self.ramses_dir, "h5_data")
        if self.success:
            self.tRelax = self.get_trelax()

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
        try:
            nml = f90nml.read(self.jobPath + "/run.sink.nml")
            t_relax = nml['POISSON_PARAMS']['trelax'] * self.unit_t / \
                          units.Myr
        except KeyError:
            t_relax = 0
        except FileNotFoundError:
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

    def overplot_time_tag(self, ax, out, timeshift=0, loc='upper left', unit="Myr",
                          **kwargs):
        """
        Overplot time tag on top-left corner
    
        Args:
            ax:
            out:
            timeshift:
    
        Returns:
    
        """
    
        pltu.overplot_time_tag(self.get_time(out) - timeshift, ax, loc=loc, unit=unit,
                          **kwargs)

    def get_out_after(self, t):
        """Get the out number right after a give time t (Myr) """
        is_first_existing = True
        for i in self.get_all_outputs():
            if self.get_time(i) > t:
                return i, is_first_existing
            is_first_existing = False
        return None, is_first_existing

    def get_gas_mass(self):
        nml = f90nml.read(self.jobPath + "/run.sink.nml")
        return nml['CLOUD_PARAMS']['mass_c']

    def get_gas_peak_density(self):
        try:
            return center.GAS_DENSITIES[self.jobid[0]]
        except KeyError:
            return -1

    def movie_sink_path(self, num):
        return "{}/movie1/sink_{:05d}.txt".format(self.jobPath, num)

    def read_movie_sink_as_particle(self, num):
        """ Read sink_xxxxx.csv in job/movie1 as particle array containing
        the following columns: m, x, y, z, vx, vy, vz, all in code units """
        
        fp = self.movie_sink_path(num)
        if not os.path.isfile(fp):
            raise FileNotFoundError
        with open(fp, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                raise NoSinkParticle
        data = np.loadtxt(fp, delimiter=',')
        if data.ndim == 1:
            data = np.array([data])
        return data[:, 1:8]

    def read_zoom_center(self):
        nml = os.path.join(self.jobPath, "run.sink.nml") 
        return read_zoom_center(nml)

    def read_zoom_radius(self, idx=-1):
        nml = os.path.join(self.jobPath, "run.sink.nml") 
        return read_zoom_radius(nml, idx)

    def get_age(self, out):
        """ return the age in Myr of all stars as an array. """

        # data = np.loadtxt("{0}/output_{1:05d}/sink_{1:05d}.csv".
        #                   format(self.jobPath, out), delimiter=',')
        data = np.loadtxt(self.get_sink_path(out), delimiter=',')
        if not len(data):
            return None
        elif len(np.shape(data)) == 1:
            data = np.array([data])
        age = data[:, -3]
        # age = age * unit['t'] / myr2s
        age = age * self.unit_t / units.Myr
        return age

    def is_sink_alive(self, out, mass_shift=1.0):
        # masses, bad = self.get_sink_mass(out)
        # if bad:
        #     return []
        try:
            masses = self.get_sink_masses(out)
        except NoSinkParticle:
            return []
        lifetime = utilities.mass_to_lifetime(mass_shift * masses)
        age = self.get_age(out)
        return age < lifetime


