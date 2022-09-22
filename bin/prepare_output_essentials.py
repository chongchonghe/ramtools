#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make the directory Job#.#/output_tmp and store into it links to all
necessary files to use yt
After this, just do the following in local computer
scp -p -r dp2:lustre/Sam/Job1.2/output_tmp/* .

@author: chongchonghe
"""

from __future__ import division, print_function
import os
import sys
import glob
import center
from tools import makedir

def link(fPath):
    """ Make sure to always work in the root directory when using this
    function to make soft link.
    """
    fName = os.path.basename(fPath)
    if not os.path.exists(fName):
        os.system("ln -s " + fPath)

def link_output(jobPath, outputID):
    tmpPath = os.path.abspath(".")
    os.chdir(jobPath)
    os.makedirs("output_tmp", exist_ok=True)
    os.chdir("output_tmp")
    outdir = "output_{:05d}".format(outputID)
    if os.path.isdir(outdir):
        return
    os.makedirs(outdir)
    os.chdir(outdir)

    out = "../output_{:05d}".format(outputID)
    link("../{}/*.info".format(out))
    link("../{}/*.csv".format(out))
    link("../{}/compilation.txt".format(out))
    link("../{}/header*.txt".format(out))
    link("../{}/hydro_file_descriptor.txt".format(out))
    link("../{}/info_0*.txt".format(out))
    link("../{}/info_rt_*.txt".format(out))
    link("../{}/namelist.txt".format(out))
    link("../{}/amr*out00001".format(out))

    os.chdir(tmpPath)

def check_and_update(job_path, is_update=True):
    """ Check and update all the soft links in job_path

    Parameters
    ----------
    job_path: string
        dir to the job path
    is_update: bool
        True: check and update soft links
        False: only check and do not update
    """
    for i in range(1, 100):
        if not os.path.isdir("{}/output_{:05d}".format(job_path, i)):
            break
        if os.path.isdir("{}/output_tmp/output_{:05d}".format(job_path, i)):
            continue
        # if INTERP goes to here, INTERP means the link to the ith output is missing
        #print("Missing {}/output_{:05d}".format(job_path, i))
        if is_update:
            print("Updating {}/output_{:05d}".format(job_path, i))
            link_output(job_path, i)
        else:
            print("Missing {}/output_{:05d}".format(job_path, i))

def check_all_jobs():
    """ *Check* all the soft links in all the jobs shown in center.py
    """
    for jobid in center.JOBIDS_ALL:
        print("Checking {} ...".format(jobid))
        check_and_update("../Job{}".format(jobid), is_update=False)

def check_and_update_all_jobs():
    """ *Check and update* all the soft links in all the jobs shown in center.py
    """
    for jobid in center.JOBIDS_ALL:
        print("Checking and updating {} ...".format(jobid))
        check_and_update("../Job{}".format(jobid))

def main(jobPath):

    os.chdir(jobPath)
    os.makedirs("output_tmp", exists_ok=True)
    os.chdir("output_tmp")
    tmpPath = os.path.abspath(".")
    for i in glob.glob("../*.nml"):
        link(i)
    for samjobs in glob.glob("../*.bash*"):
        link(samjobs)
    for out in glob.glob("../output_0*"):
        outbase = os.path.basename(out)
        if os.path.islink(out):
            if os.path.isdir(outbase):
                os.system('rm -rf {}'.format(outbase))
            if not os.path.isfile(outbase + '.link'):
                os.system("echo '{}' > {}.link".format(os.readlink(out), outbase))
            continue
        if os.path.isdir(outbase):
            continue
        os.makedirs(outbase)
        os.chdir(outbase)
        link("../{}/*.info".format(out))
        link("../{}/*.csv".format(out))
        link("../{}/compilation.txt".format(out))
        link("../{}/header*.txt".format(out))
        link("../{}/hydro_file_descriptor.txt".format(out))
        link("../{}/info_0*.txt".format(out))
        link("../{}/info_rt_*.txt".format(out))
        link("../{}/namelist.txt".format(out))
        link("../{}/amr*out00001".format(out))
        os.chdir(tmpPath)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: {} path".format(sys.argv[0]))
    #check_and_update_all_jobs()
    main(sys.argv[1])
