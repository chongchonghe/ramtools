import os
import numpy as np
import yt
from .yt_field_descrs import FIELDS

def mass_to_temp(mass):
    """ input: mass (M_sun); output: T (K) """

    radius = 1.06 * mass**0.945 if mass < 1.66 else 1.33 * mass**0.555
    if mass > 20.0:
        lumino = 81 * mass**2.14
    elif mass > 2.0:
        lumino = 1.78 * mass**3.5
    else:
        lumino = 0.75 * mass**4.8
    return 5777 * (lumino / radius**2)**(0.25)

TLIM = [mass_to_temp(10 ** -1), mass_to_temp(10 ** 3)]

def my_yt_load(job_path, out):
    """Quick load a job to yt with the correct FIELDS

    Args:
        job_path (str): path to a job directory
        out (int): the output you want to load

    """

    return yt.load("{}/output_{:05d}/info_{:05d}.txt".format(
        job_path, out, out), fields=FIELDS)

def read_quant_from_ramses_info(info, quant):
    with open(info, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line == '':
            continue
        words = line.split('=')
        var = words[0].strip()
        if ' ' in var:
            break
        if var == quant:
            return float(words[1].strip())
    return None


def get_times_from_movie1(movie_dir):
    """ Return a dictionary: {100: 0.012, 101: 0.013, ...}
    """

    sink_fn_fmt = movie_dir + "/sink_{:05d}.txt"
    info_fn_fmt = movie_dir + "/info_{:05d}.txt"
    data = {}
    for i in range(5000):
        fn = sink_fn_fmt.format(i)
        if not os.path.isfile(fn):
            continue
        with open(fn, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                continue
        info = info_fn_fmt.format(i)
        t = read_quant_from_ramses_info(info, 'time')
        data[i] = t
    return data


def get_sink_info_from_movie1(movie_dir, sinkid):
    """ Return a list of time and a list of sink masses

    Args
        sinkid: interger, starting from 0

    Return
    ------
    Return a dictionary containing the following keys 'out', 't', 'm', 'x',
    'y', 'z', 'vx', 'vy', 'vz'. All in code units.

    """

    sink_fn_fmt = movie_dir + "/sink_{:05d}.txt"
    info_fn_fmt = movie_dir + "/info_{:05d}.txt"
    data = {}
    elems = ['out', 't', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    for elem in elems:
        data[elem] = []
    for i in range(5000):
        fn = sink_fn_fmt.format(i)
        if not os.path.isfile(fn):
            continue
        with open(fn, 'r') as f:
            if not os.fstat(f.fileno()).st_size:
                continue
            num_lines = len(f.readlines())
        if num_lines - 1 < sinkid:
            continue
        pars_new = np.loadtxt(fn, delimiter=',')
        if pars_new.ndim == 1:
            pars_new = np.array([pars_new])
        parnew = pars_new[sinkid, :]
        data['out'].append(i)
        data['m'].append(parnew[1])
        data['x'].append(parnew[2])
        data['y'].append(parnew[3])
        data['z'].append(parnew[4])
        data['vx'].append(parnew[5])
        data['vy'].append(parnew[6])
        data['vz'].append(parnew[7])

        fn = info_fn_fmt.format(i)
        t = read_quant_from_ramses_info(fn, 'time')
        data['t'].append(t)
    return data
