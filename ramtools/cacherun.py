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
from collections import OrderedDict
import yt
from .utilities import my_yt_load
try:
    yt.set_log_level(50)
except:
    pass

CACHEPATH = None
GLOBALFORCEREPLACECACHE = False

class CacheFile(object):

    def __init__(self, snapshot, algorithm):
        self._snapshot = snapshot  # the absolute path to an output folder
        self._algorithm = algorithm  # the function to run on this snapshot
        self._folder = snapshot.CachePath()
    
    def Save(self, data):
        # Save algorithm settings to a text file (for reference)
        pikfile = open(self._Filename("info"),"wb")
        pik.dump(str(self._algorithm),pikfile)
        pikfile.close()
        # Save the output data to a binary file
        pikfile = open(self._Filename("data"),"wb")
        pik.dump(data,pikfile)
        pikfile.close()
        
    def Load(self):
        # Load the (binary) data file
        if self.Exists():
            pikfile = open(self._Filename("data"),"rb")
            print("Loading from cache: snap", self._snapshot._outnum,\
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
        return self._folder+self._algorithm.CacheFilename(ext)

class HashRun(object):
    '''
    A wrapper class around a function that enables us to run functions and store their outputs for later use
    '''

    def __init__(self, function):
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
        self._force_replace_cache = False

    def ForceReplaceCache(self,on):
        self._force_replace_cache = on

    def FunctionName(self):
        '''
        Return the function's name
        '''
        return self._function.__name__
        
    def Run(self, snap, *args, **kwargs):
        '''
        Run for a single snapshot
        '''
        global GLOBALFORCEREPLACECACHE
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(snap, self)
        # Check to see if a cached dataset exists
        if cache.Exists() and not self._force_replace_cache and not GLOBALFORCEREPLACECACHE:
            try:
                output = cache.Load()
            except EOFError:
                # If the cache is corrupted, rerun anyway
                output = self._RunAlgorithm(snap)
                cache.Save(output)
        else:
            output = self._RunAlgorithm(snap)
            cache.Save(output)
        return output

    def Cached(self, snap, *args, **kwargs):
        '''
        Check that the cache for this data exists
        '''
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(snap, self)
        # Check to see if a cached dataset exists
        return cache.Exists()
        
    def _RunAlgorithm(self, snapshot):
        '''
        Unpack the algorithm and call the native python code
        '''
        ds = snapshot._ds
        output = self._function(ds, *self._args, **self._kwargs)
        return output
    
    def __str__(self):
        '''
        Parse this algorithm as a string (for writing to text file)
        '''
        out = ""
        out += "Function: "+self._function.__name__
        out += "\n"
        out += "Arguments: "+str(self._args)
        out += "\n"
        out += "Keyword Arguments: "+str(self._kwargs)
        out += "\n"
        return out
    
    def __call__(self, snapshot, *args, **kwargs):
        '''
        Allows the algorithm to be called like a function
        '''
        return self.Run(snapshot, *args, **kwargs)
    
    def CacheFilename(self, ext="data"):
        '''
        Return the name of this algorithm's cache file
        Uses a hash function to hash the arguments of the function
        ext - File extension (used for writing multiple cache files; default is ".data")
        '''
        objName = self._function.__name__
        hash = hashlib.sha1(str(self._args).encode('utf-8')).hexdigest()
        hash += hashlib.sha1(str(self._kwargs).encode('utf-8')).hexdigest()
        filepre = objName+hash
        return filepre+"."+ext

class RamsesSnapshot(object):

    def __init__(self,folder,outnum,name=None):
        self._ds = my_yt_load(folder, outnum) # yt.load with correct FIELDS
        self._folder = folder
        self._outnum = outnum
        self._name = name

    # def _Setup(self):
    #     if self._ds is None:
    #         # self._ds = pymses.RamsesOutput(self._folder,self._outnum)
    #         self._ds = load_ds(self._folder, self._outnum)
    #         # # Hack to allow nested Algorithm calls
    #         # self._ds.hamusnap = self

    # def RawData(self):
    #     self._Setup()
    #     return self._ds

    # def Name(self):
    #     return self._name

    def CachePath(self):
        #if self._name is None:
        #    pathname = self._ds.output_repos
        #else:
        #    pathname = self._name
        #pathname = pathname.replace("/","__").replace("~","TILDE")
        #path = "./cache/"+pathname+"/output_"+str(self.OutputNumber()).zfill(5)+"/"
        if CACHEPATH is not None:
            pathname = CACHEPATH+"/"+self._folder
        else:
            pathname = self._folder
        # Put the cache files in the simulation folder itself
        path = f"{pathname}/cache/yt{yt.__version__}/output_{self._outnum:05d}/"
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except:
                pass # Probably fine. Probably.
        return path
