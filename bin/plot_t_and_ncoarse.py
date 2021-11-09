#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" timeOutFile.py
Read Job*.bash.o* file and get_t_and_ncoarse

Usage:
>>> get_t_and_ncoarse(nml_fielname)
"""

import sys, os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
from glob import glob

def get_t_and_ncoarse(fn):
    """Given output file name fn, return the list of n_coarse and
    dt_coarse
    """

    togrep = "Time elapsed since last"
    try:
        ret = subprocess.check_output("grep '{}' {}".format(togrep, fn), shell=True)
    except subprocess.CalledProcessError:
        print("Called Error")
        return [float('nan')], [float('nan')]
    ret = ret.decode()
    words = ret.split()
    switch = 0
    ts = []
    for word in words:
        if word == 'step:':
            switch = 1
            continue
        if switch:
            ts.append(float(word))
            switch = 0
        if len(word) > 5:
            if word[:5] == 'step:':
                ts.append(float(word[5:]))

    try:
        ret = subprocess.check_output("grep -m 1 'nstep_coarse' {}".format(fn),
                                      shell=True)
        nc = ret.decode().split()[5]
        nc = int(nc)
    except subprocess.CalledProcessError:
        nc = 0
    return list(np.arange(nc, nc + len(ts), dtype=int)), ts


def get_info_arg(fn, arg):
    with open(fn, 'r') as f:
        for line in f.readlines():
            words = line.split()
            if len(words) >= 2:
                if words[0] in [arg, arg + '=']:
                    return float(line.split()[-1])
    print('Failed to read arg in', fn)
    return None


def plot_dt_vs_coarse(jobdir, ax=None):

    if ax is None:
        plt.figure()
        ax = plt.gca()

    ncs = []
    tss = []
    for f in sorted(glob(jobdir + '/slurm-*.out')):
        nc, ts = get_t_and_ncoarse(f)
        lb = os.path.basename(f)
        ax.plot(nc, np.array(ts) / 60, '.', label=lb)
        ncs += nc
        tss += ts
        print('In {} found n_coarse_i = {}'.format(f, nc[0]))

    # figure out how long is a coarse step
    outs = sorted(glob('jobdir' + '/output_0*/info_0*.txt'))
    if len(outs) == 0:
        return
    ts = []
    ns = []
    for out in outs:
        ts.append( get_info_arg(out, 'time') )
        ns.append( get_info_arg(out, 'nstep_coarse') )
    ts = np.array(ts)
    ns = np.array(ns)
    dts = (ts[1:] - ts[:-1]) / (ns[1:] - ns[:-1])
    mean = np.mean(dts)
    std = np.std(dts)
    myr_to_s = 3.154e13
    try:
        unit_t = get_info_arg(outs[0], 'unit_t')
    except IndexError:
        print(outs)
        return
    t_c_to_kyr = unit_t / myr_to_s * 1e3

    # place a text box in upper left in axes coords
    textstr = '\n'.join([
        't_coarse = {:.2e} = {:.1f} kyr'.format(mean, mean * t_c_to_kyr),
        'std = {:.2e} = {:.1f} kyr'.format(std, std * t_c_to_kyr)])
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize='large',
            va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel('n_coarse')
    ax.set_ylabel('Time elapsed (m) in each coarse step')
    tt = os.path.basename(os.path.abspath('.'))
    ax.set_title(tt)
    ax.legend(loc='center left')
    # plt.savefig('coarse_time.png', dpi=200)

def main():

    _, ax = plt.subplots(figsize=[12, 5])
    plot_dt_vs_coarse('.', ax)
    plt.title(os.path.basename(os.path.abspath('.')))
    plt.tight_layout()
    plt.savefig('t_and_ncoarse.png')

if __name__ == "__main__":

    main()
    if sys.platform == "darwin":
        os.system('open t_and_ncoarse.png')
