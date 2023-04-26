"""
Run YT functions with cache support so that when a function is run a second
time, a cache from the first run is used to return the results immediately.

Adapted from Sam Geen's HamuLite.py which describes as 'MODIFIED FROM
Algorithm.py IN HAMU TO WORK INDEPENDENTLY OF THE HAMU FRAMEWORK'

Original script:
https://github.com/samgeen/mcrtscripts/blob/master/scripts/HamuLite.py
"""

import sys, os
if sys.version_info[0] < 3:
    import cPickle as pik
else:
    import pickle as pik
import hashlib
import warnings
from collections import OrderedDict
import yt
from .utilities import my_yt_load
try:
    yt.set_log_level(50)
except:
    pass

CACHEPATH = None
GLOBALFORCEREPLACECACHE = False
ISLOG = True

def turn_on_global_force_replace_cache():
    global GLOBALFORCEREPLACECACHE
    GLOBALFORCEREPLACECACHE = True

def turn_off_logging():
    global ISLOG
    ISLOG = False

def get_jobdir_and_out_from_ds(ds):
    """Given ds, which has ds.directory = '/path/Job1/output_00002', return
    /path/Job1 and output_00002 """
    outdir = os.path.abspath(ds.directory)
    jobdir = os.path.dirname(outdir)
    return jobdir, os.path.basename(outdir)


class CacheFile(object):

    def __init__(self, ds, algorithm):
        self._snapshot = ds  # the absolute path to an output folder
        self._algorithm = algorithm  # the function to run on this snapshot
        jobdir, out = get_jobdir_and_out_from_ds(ds)
        self._folder = f"{jobdir}/cache/cache/yt{yt.__version__}/{out}"
        if not os.path.exists(self._folder):
            os.makedirs(self._folder)

    def Save(self, data):
        # Save algorithm settings to a text file (for reference)
        pikfile = open(self._Filename("info"),"wb")
        pik.dump(str(self._algorithm),pikfile)
        pikfile.close()
        # Save the output data to a binary file
        fo = self._Filename("data")
        pikfile = open(fo, "wb")
        pik.dump(data,pikfile)
        pikfile.close()
        if ISLOG:
            print(f"Saving cache file {fo}.")
        
    def Load(self):
        # Load the (binary) data file
        if self.Exists():
            pikfile = open(self._Filename("data"),"rb")
            if ISLOG:
                print("Loading from cache: ds", self._folder,\
                    self._algorithm.FunctionName(), "...")
            try:
                output = pik.load(pikfile)
            except:
                print("Hamu: Error reading cache file", pikfile)
                raise ValueError
            pikfile.close()
            return output
        else:
            return None
        
    def Exists(self):
        '''
        Does the cache file exist?
        '''
        return os.path.exists(self._Filename())
        
    def _Filename(self,ext="data"):
        '''
        Cache file's filename
        '''
        return os.path.join(self._folder, self._algorithm.CacheFilename(ext))


class CacheRun(object):
    '''
    A wrapper class around a function that enables us to run functions and store their outputs for later use

    Usage (new):
    >>> def get_max_pos(ds, field):
    >>>     return ds.all_data().argmax(field)
    >>> CacheRun(get_max_pos)(ds, ('gas', 'density'))

    Usage (old):
    >>> ds = RamsesSnapshot("tests/Job2", 40)
    >>> def get_max_pos(ds, field):
    >>>     return ds.all_data().argmax(field)
    >>> CacheRun(get_max_pos)(ds, ('gas', 'density'))
    '''

    def __init__(self, function, force_replace_cache=False):
        '''
        # def __init__(self, function, *args, **kwargs):
        Constructor
        function: A Python function object to call that accepts snapshot._ds
        as its first argument, and arg/kwarg after
        arg, kwarg: arguments accepted by function (see, e.g., http://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/)
        '''
        self._function = function
        # self._args = args
        # self._kwargs = kwargs
        self._force_replace_cache = force_replace_cache

    def ForceReplaceCache(self, on=True):
        self._force_replace_cache = on

    def FunctionName(self):
        '''
        Return the function's name
        '''
        return self._function.__name__

    def Run(self, ds, *args, **kwargs):
        '''
        Run for a single snapshot
        '''
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(ds, self)
        # Check to see if a cached dataset exists
        if cache.Exists() and not self._force_replace_cache and not GLOBALFORCEREPLACECACHE:
            try:
                output = cache.Load()
            except EOFError:
                # If the cache is corrupted, rerun anyway
                warnings.warn(f"{cache._Filename()} failed to load. Will run the "
                              f"function explicitly")
                output = self._RunAlgorithm(ds)
                cache.Save(output)
        else:
            output = self._RunAlgorithm(ds)
            cache.Save(output)
        return output

    # def Cached(self, snap, *args, **kwargs):
    #     '''
    #     Check that the cache for this data exists
    #     '''
    #     self._args = args
    #     self._kwargs = kwargs
    #     # First, get the cache filename and compare against existing files
    #     cache = CacheFile(snap, self)
    #     # Check to see if a cached dataset exists
    #     return cache.Exists()

    def _RunAlgorithm(self, snapshot):
        '''
        Unpack the algorithm and call the native python code
        '''
        ds = snapshot
        output = self._function(ds, *self._args, **self._kwargs)
        return output

    def __str__(self):
        '''
        Parse this algorithm as a string (for writing to text file)
        '''
        out = ""
        out += "Function: " + self._function.__name__
        out += "\n"
        out += "Arguments: " + str(self._args)
        out += "\n"
        out += "Keyword Arguments: " + str(self._kwargs)
        out += "\n"
        return out

    def __call__(self, ds, *args, **kwargs):
        '''
        Allows the algorithm to be called like a function
        '''
        return self.Run(ds, *args, **kwargs)

    def CacheFilename(self, ext="data"):
        '''
        Return the name of this algorithm's cache file
        Uses a hash function to hash the arguments of the function
        ext - File extension (used for writing multiple cache files; default is ".data")
        '''
        objName = self._function.__name__
        hash = hashlib.sha1(str(self._args).encode('utf-8')).hexdigest()
        # hash += hashlib.sha1(str(self._kwargs).encode('utf-8')).hexdigest()
        cleaned_kwargs = {ite:(self._kwargs[ite] if not callable(
            self._kwargs[ite]) else self._kwargs[ite].__name__) for ite in
                          self._kwargs.keys()}
        hash += hashlib.sha1(str(cleaned_kwargs).encode('utf-8')).hexdigest()
        filepre = objName + '-' + hash
        return filepre + "." + ext
