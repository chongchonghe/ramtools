import yt
from .yt_field_descrs import FIELDS

def my_yt_load(job_path, out):
    """Quick load a job to yt with the correct FIELDS

    Args:
        job_path (str): path to a job directory
        out (int): the output you want to load

    """

    return yt.load("{}/output_{:05d}/info_{:05d}.txt".format(
        job_path, out, out), fields=FIELDS)

