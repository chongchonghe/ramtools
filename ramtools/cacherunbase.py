"""
A general protocal for storing function result in cache files

What is does: The first time you run CacheRun, the result is stored in a file using pickle and is named
after a unique hash string based on
(1) the path of the Python file you are running,
(2) the function name, and
(3) the parameters passed to the function.
The second time you run CacheRun with the same hash, the previously stored data is loaded and returned.

Example:
>>> def func(x, y=1):
>>>     return x + y
>>> CacheRun(func, os.path.abspath(__file__))(10, y=2)
>>> CacheRun(func, os.path.abspath(__file__))(10, y=2)
"""

import sys, os
if sys.version_info[0] < 3:
    import cPickle as pik
else:
    import pickle as pik
import hashlib
import warnings

CACHEPATH = None
GLOBALFORCEREPLACECACHE = False


def turn_on_global_force_replace_cache():
    global GLOBALFORCEREPLACECACHE
    GLOBALFORCEREPLACECACHE = True


class CacheFile(object):

    def __init__(self, folder, algorithm):
        self._algorithm = algorithm  # the function to run on this snapshot
        self._folder = os.path.join(folder, 'cache')
        if not os.path.exists(self._folder):
            os.makedirs(self._folder)

    def Save(self, data):
        # Save algorithm settings to a text file (for reference)
        pikfile = open(self._Filename("info"), "wb")
        pik.dump(str(self._algorithm), pikfile)
        pikfile.close()
        # Save the output data to a binary file
        fo = self._Filename("data")
        pikfile = open(fo, "wb")
        pik.dump(data, pikfile)
        pikfile.close()
        print(f"Saving cache file {fo}.")

    def Load(self):
        # Load the (binary) data file
        if self.Exists():
            pikfile = open(self._Filename("data"), "rb")
            print("Loading from cache: ", self._folder, \
                  self._algorithm._function.__name__, "...")
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

    def _Filename(self, ext="data"):
        '''
        Cache file's filename
        '''
        return os.path.join(self._folder, self._algorithm.CacheFilename(ext))


class CacheRun():

    def __init__(self, function, filename, force_replace_cache=False, flag=''):
        """

        Args:
            function: the function to run
            filename: == os.path.abspath(__file__)
            force_replace_cache: toggle force replace cache
            flag: a unique flag to tell apart functions with same names
                inside different functions, or different versions of the
                function, although in the latter case a force_replace_cache is
                recommended.

        Examples:
            >>> def func(x, y=1):
            >>>     return x + y
            >>> CacheRun(func, os.path.abspath(__file__))(10, y=2)
        """

        self._function = function
        self._filename = filename
        self._flag = flag
        self._force_replace_cache = force_replace_cache

    def ForceReplaceCache(self, on=True):
        self._force_replace_cache = on

    def Run(self, flag, *args, **kwargs):
        '''
        '''

        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(self)
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

    def CacheFilename(self, ext="data"):
        '''
        Return the name of this algorithm's cache file
        Uses a hash function to hash the arguments of the function
        ext - File extension (used for writing multiple cache files; default is ".data")
        '''
        filename = os.path.basename(self._filename)
        objName = self._function.__name__
        hash = hashlib.sha1(str(self._args).encode('utf-8')).hexdigest()
        # hash += hashlib.sha1(str(self._kwargs).encode('utf-8')).hexdigest()
        cleaned_kwargs = {ite:(self._kwargs[ite] if not callable(
            self._kwargs[ite]) else self._kwargs[ite].__name__) for ite in
                          self._kwargs.keys()}
        cleaned_kwargs['flag'] = self._flag
        hash += hashlib.sha1(str(cleaned_kwargs).encode('utf-8')).hexdigest()
        filepre = filename + '-' + objName + '-' + hash
        return filepre + "." + ext

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

    def _RunAlgorithm(self):
        '''
        Unpack the algorithm and call the native python code
        '''
        output = self._function(*self._args, **self._kwargs)
        return output

    def Run(self, *args, **kwargs):
        '''
        Run for a single snapshot
        '''
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cachefolder = os.path.dirname(self._filename)
        cache = CacheFile(cachefolder, self)
        # Check to see if a cached dataset exists
        if cache.Exists() and not self._force_replace_cache and not GLOBALFORCEREPLACECACHE:
            try:
                output = cache.Load()
            except EOFError:
                # If the cache is corrupted, rerun anyway
                warnings.warn(f"{cache._Filename()} failed to load. Will run the "
                              f"function explicitly")
                output = self._RunAlgorithm()
                cache.Save(output)
        else:
            output = self._RunAlgorithm()
            cache.Save(output)
        return output

    def __call__(self, *args, **kwargs):
        '''
        Allows the algorithm to be called like a function
        '''
        return self.Run(*args, **kwargs)

