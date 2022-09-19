
import os
import warnings
import json
import numpy as np
from matplotlib.pyplot import cm
import f90nml
import yt
try:
    yt.set_log_level(50)
except:
    pass

from . import utilities as util
from . import units


class NoSinkParticle(Exception):
    """exception: no sink particle found in an output or a .csv file."""
    pass


def set_RAM_DIR(ram_dir):
    """Set the global variable RAM_DIR. Check Ramses.__init__ for its usage."""

    global RAM_DIR
    RAM_DIR = ram_dir


def to_boxlen(quant, ds):
    if type(quant) == tuple:
        return float(yt.YTQuantity(*quant) / ds.length_unit)
    else:
        return quant


class RamsesBase():
    """This is the core of ramtools. Most of the times, you want to start ramtools
by declaring this class.

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

    def __init__(self, jobdir, fields=None):
        """
        Parameters
        ----------
        jobdir: str 
            absolute/relative path to the job 
        dirname: str (default None)
            the basename of the job directory. You need to do ramtools.set_RAM_DIR("path/to/a/folder/of/RAMSES/jobs")
        fields: dict
            A dictionary of RAMSES FIELDS

        Returns
        -------
        0 if successfully loaded a ds; 1 otherwise

        """
        self.jobPath = jobdir
        self.fields = fields
        self.ds_container = {}
        self.ds1 = None
        self.get_ds()
        if self.get_units():
            self.success = False
        else:
            self.success = True

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
        if self.ds1 is None:
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
        # self.unit_B = util.get_unit_B(self.ds1)
        self.unit_B = util.get_unit_B_new(self.ds1)
        # self.n_colden_H = np.double(self.unitDen2Hcc * self.unit_l)
        return 0

    def get_info_path(self, out):
        """Return the path to info_out.txt"""
        return "{0}/output_{1:05d}/info_{1:05d}.txt".format(
            self.jobPath, out)

    def exist(self, out):
        return os.path.exists(self.get_info_path(out))

    def load_ds(self, out):
        """Return a yt.load instance of the frame `out`"""
        if out not in self.ds_container:
            if self.fields is None:
                self.ds_container[out] = yt.load(self.get_info_path(out))
            else:
                self.ds_container[out] = yt.load(self.get_info_path(out), fields=self.fields)
        return self.ds_container[out]

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
        for i in range(1, 200):
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
        for i in range(200, 0, -1):
            fname = "{0}/output_{1:05d}/info_{1:05d}.txt".format(
                self.jobPath, i)
            if os.path.isfile(fname):
                return i
        return 0
