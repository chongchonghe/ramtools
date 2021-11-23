"""
A general wrapper to make some fo the common yt functions more efficient by reusing intermediate data
"""

import os
import yt
from typing import Dict, Any
import hashlib
import json

try:
    yt.set_log_level(40)
except:
    pass

DATA_DIR = '.'

def dict_hash(dictionary: Dict[str, Any]) -> str:
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # Sort arguments so {'a': 1, 'b': 2} is the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()

def set_data_dir(datadir):
    global DATA_DIR
    DATA_DIR = datadir

def to_tuple(x):
    if isinstance(x, list):
        return tuple(x)
    return x

def ProfilePlot(ds, axis, fields, center='c', width=None,
                axes_unit=None, weight_field=None, max_level=None,
                origin='center-window', tag=None, **kwargs):

    if isinstance(axis, int):
        axis = ['x', 'y', 'z'][axis]
    hash_dict = dict(
        ds = str(ds.__repr__()),
        axis = axis,
        fields = fields,
        center = to_tuple(center),
        width = width, 
        weight_field = weight_field,
        max_level = max_level,
        kind = "prj_plot",
        tag = tag,  # a free tag, used to distinguish different versions 
    )
    # print(hash_dict)
    hashstr = dict_hash(hash_dict)
    h5fn = os.path.join(DATA_DIR, hashstr + ".h5")
    if not os.path.exists(h5fn):
        if not os.path.isdir(DATA_DIR):
            os.makedirs(DATA_DIR)
        p = ds.proj(field=fields, axis=axis, center=center,
                    weight_field=weight_field,
                    field_parameters={'width': width},
                    )
        p.save_as_dataset(h5fn)
    else:
        print("Loading", h5fn)
    data = yt.load(h5fn)
    p = yt.ProjectionPlot(data, axis, fields, center=center,
                          width=width, weight_field=weight_field,
                          axes_unit=axes_unit, max_level=max_level, 
                          **kwargs)
    return p

def SlicePlot():
    pass