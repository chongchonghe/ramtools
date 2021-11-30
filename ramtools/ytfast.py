"""A general wrapper to make some fo the common yt functions more
efficient by reusing intermediate data

General-purpose yt wrappers. Apply to any data that yt supports.

"""

import os
import numpy as np
import yt
from matplotlib import colors
import matplotlib.pyplot as plt
from typing import Dict, Any
import hashlib
import json
from time import time

try:
    import logging
    logger = logging.basicConfig(level=logging.INFO)
    Logger = logging.getLogger(__name__)
    Logger.setLevel(logging.ERROR)
    ISLOG = True
except ModuleNotFoundError:
    ISLOG = False
    pass

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

def ProjectionPlot(ds, axis, fields, center='c', width=None,
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

def PhasePlot(data_source, x_field, y_field, z_fields,
              weight_field=None, x_bins=128, y_bins=128,
              accumulation=False, fractional=False, fontsize=18,
              figure_size=8.0, shading='nearest', extrema=None,
              units=None, zlims=None, force_redo=False,
              define_field=None, is_cb=True, cmap="viridis",
              cb_label='', f=None, ax=None, ret='imshow',):
    """

    Args:
        # args passed to yt.PhasePlot
        data_source (YTSelectionContainer Object)
        x_field (tuple)
        y_field (tuple)
        z_fields (list of tuples)
        weight_field (tuple)
        # args passed to yt.PhasePlot
        zlims (list_like): lims of z_fields
        force_redo (bool): toggle always make h5 data
        define_field (func): function to define a new field
        is_cb (bool): toggle show colorbar (default: True)
        cb_label (str): colorbar label (default: '')
        f (plt.figure): figure to plot on (default: None)
        ax (plt.axis): axis to plot on (default: None)
        ret (str): what to return (default: 'imshow'). One of ['imshow', 'data']

    """

    # try:
    #     data_source.ds
    # except AttributeError:
    #     Warning("data_source has to be YTSelectionContainer Object, e.g. all_data, sphere")
    #     return
    # logs_str = {str(key): value for key, value in logs.items()}

    t1 = time()
    if not isinstance(x_field, tuple):
        x_field = ('gas', x_field)
    if not isinstance(y_field, tuple):
        y_field = ('gas', y_field)
    if not isinstance(z_fields[0], tuple):
        z_fields[0] = ('gas', z_fields[0])
    if extrema is None:
        extrema = {x_field: [None, None], y_field: [None, None]}
    extrema_str = {str(key): value for key, value in extrema.items()}
    # if units is not None:
    #     units_str = {str(key): value for key, value in units.items()}
    # else:
    #     units_str = units

    t1 = time()
    hash_params = dict(
        ds = str(data_source.ds.__repr__()),
        data_source = str(data_source.__repr__()),
        x_field = x_field, y_field = y_field, z_fields = z_fields,
        weight_field = weight_field, x_bins=x_bins, y_bins=y_bins,
        extrema = extrema_str, accumulation=accumulation,
    )
    hashstr = dict_hash(hash_params)
    h5fn = os.path.join(DATA_DIR, hashstr + ".h5")
    jsonfn = os.path.join(DATA_DIR, hashstr + ".json")
    if ISLOG: Logger.info(f"Log 1, dt = {time() - t1}")
    if force_redo or (not os.path.exists(h5fn)):
        if define_field is not None:
            define_field(data_source.ds)
        assert x_bins == y_bins
        p = yt.create_profile(data_source, [x_field, y_field],
                              z_fields, weight_field=weight_field,
                              extrema=extrema, n_bins=x_bins,
                              )
        p.save_as_dataset(h5fn)
        with open(jsonfn, 'w') as fi:
            json.dump(hash_params, fi, indent=2)
    if ISLOG: Logger.info(f"Log 2, dt = {time() - t1}")
    prof = yt.load(h5fn)
    if ret == 'data':
        return prof
    if ISLOG: Logger.info(f"Log 3, dt = {time() - t1}")
    # print(f"time: {time() - t1}")
    dat = prof.data[z_fields[0]].T
    extents = [prof.data[x_field].min(), prof.data[x_field].max(),
                prof.data[y_field].min(), prof.data[y_field].max()]
    extent_log = np.log10(extents)
    if zlims is None:
        thenorm = colors.LogNorm()
    else:
        thenorm = colors.LogNorm(vmin=zlims[0], vmax=zlims[1])
    if f is None and ax is None:
        f, ax = plt.subplots()
    if ISLOG: Logger.info(f"Log 4, dt = {time() - t1}")
    p = ax.imshow(dat,
                  norm=thenorm,
                  extent=extent_log,
                  aspect="auto", origin="lower",
                  cmap=cmap,
                  )
    if ISLOG: Logger.info(f"Log 5, dt = {time() - t1}")
    if is_cb and f is not None:
        cb = f.colorbar(p)
        cb.set_label(cb_label)
    str1 = x_field[1] if isinstance(x_field, tuple) else x_field
    str2 = y_field[1] if isinstance(y_field, tuple) else y_field
    ax.set(xlabel="log " + str1, ylabel=r"log " + str2)
    del dat
    if ISLOG: Logger.info(f"Log 6, dt = {time() - t1}")
    if is_cb:
        return f, ax, p, cb
    else:
        return f, ax, p

def ProfilePlot(data_source, x_field, y_fields, weight_field=('gas', 'mass'),
                n_bins=64, accumulation=False, fractional=False, label=None,
                plot_spec=None, x_log=True, y_log=True, xlims=[None, None],
                force_redo=False, define_field=None, f=None, ax=None,
                ret='plot', mpl_kwargs={}):

    hash_params = dict(
        ds = str(data_source.ds.__repr__()),
        data_source = str(data_source.__repr__()),
        x_field = x_field, y_fields = y_fields,
        weight_field = weight_field,
    )
    if n_bins != 64:
        hash_params['n_bins'] = n_bins
    if accumulation:
        hash_params['accumulation'] = accumulation
    if fractional:
        hash_params['fractional'] = fractional
    # if not x_log:
    #     hash_params['x_log'] = x_log
    # if not y_log:
    #     hash_params['y_log'] = y_log
    if not xlims == [None, None]:
        hash_params['xlims'] = xlims
    hashstr = dict_hash(hash_params)
    h5fn = os.path.join(DATA_DIR, '1dprofile-' + hashstr + ".h5")
    jsonfn = os.path.join(DATA_DIR, '1dprofile-' + hashstr + ".json")
    if force_redo or (not os.path.exists(h5fn)):
        print(f"Making {h5fn}")
        if define_field is not None:
            define_field(data_source.ds)
        p = yt.create_profile(data_source, x_field, y_fields,
                              weight_field=weight_field, n_bins=n_bins,
                              extrema={x_field: xlims}, accumulation=accumulation,
                              fractional=fractional,
                              )
        p.save_as_dataset(h5fn)
        with open(jsonfn, 'w') as fi:
            json.dump(hash_params, fi, indent=2)
    prof = yt.load(h5fn)
    x = prof.data[x_field]
    y = prof.data[y_fields[0]]
    if ret == 'data':
        return x, y
    if f is None and ax is None:
        f, ax = plt.subplots()
    ax.plot(x, y, **mpl_kwargs)
    if x_log: ax.set_xscale('log')
    if y_log: ax.set_yscale('log')
    str1 = x_field[1] if isinstance(x_field, tuple) else x_field
    str2 = y_fields[0][1] if isinstance(y_fields[0], tuple) else y_fields[0]
    if x_log: str1 = "log " + str1
    if y_log: str2 = "log " + str2
    ax.set(xlabel=str1, ylabel=str2)
    return f, ax
