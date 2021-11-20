# -*- coding: utf-8 -*-
"""
Some utilities to make figures.

Attributes:
    ...
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

from .ramses import Ramses

def dict_hash(dictionary: Dict[str, Any]) -> str:
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # Sort arguments so {'a': 1, 'b': 2} is the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


def den_setup(p, zlim=None, time_offset=None, mu=1.4, unit='number_density',
              weight_field=None):
    """
    Args:
        p (YT plot): 
        zlim (tuple): limits of the field
        mu (float): the mean particle weight. The mass density rho = mu * n
        quant (str): 'volume' or 'column' for volumetric and column density
    """

    if unit == 'number_density':
        if weight_field == "density":
            p.set_unit(('gas', 'density'), 'cm**-3',
                    equivalency='number_density',
                    equivalency_kwargs={'mu': mu})
    if zlim is not None:
        p.set_zlim(('gas', 'density'), zlim[0], zlim[1])
    # This will change the pixel size of the figure. Be careful!
    p.set_figure_size(6)
    p.set_cmap(('gas', 'density'), 'inferno')
    p.annotate_timestamp(time_format='{time:.3f} {units}',
                         time_offset=time_offset)


def T_setup(p, zlim=None, time_offset=None):
    if zlim is not None:
        p.set_zlim('temperature', zlim[0], zlim[1])
    p.set_figure_size(6)
    p.set_cmap('temperature', 'gist_heat')
    p.annotate_timestamp(time_format='{time:.3f} {units}',
                         time_offset=time_offset)


class RamPlot(Ramses):

    def to_boxlen(self, quant):
        if isinstance(quant, tuple):
            return yt.YTQuantity(*quant).to('cm').value / self.unit_l
        else:
            return quant

    def plot_prj(self,
                 out,
                 ds=None,
                 center=(0.5, 0.5, 0.5),
                 width=0.5,
                 axis='x',
                 field="density",
                 weight_field="density",
                 kind="projection",
                 normal=None,
                 north=None,
                 tag=0,
                 ):
        """Plot projection or slice plot. When plotting projection plot, an 
        intermediate data is stored for fast recreate of the figure.

        Args:
            out (int): the frame number
            ds (YT ds): optional. If not specified, will load a ds using self.load_ds
            center (tuple): the center as a tuple in boxlen unit. Default: (0.5, 0.5, 0.5)
            width (tuple or float): the width of the area to plot. e.g. (0.1 'pc'); 0.02. Default: 0.5
            axis (int or str): the axis of the line of sight. Default: 'x'
            field (tuple or str): the field to plot. Default: 'density'
            weight_field (tuple or str): the field to weight with. Default: 'density'
            kind (str): the figure type, either 'projection' (default) or 'slice'.
            normal (tuple or list): 
        
        Returns:
            p (YT plot instance): you can save it to file by `p.save('filename.png')`

        """

        width = self.to_boxlen(width)
        if isinstance(axis, int):
            axis = ['x', 'y', 'z'][axis]
        if isinstance(field, str):
            if field in ["density", "temperature", "pressure"]:
                field = ("gas", field)
        if isinstance(center, list):
            center = tuple(center)
        if isinstance(normal, list):
            normal = tuple(normal)
        if isinstance(north, list):
            north = tuple(north)

        params = dict(
            jobdir = os.path.abspath(self.jobPath),
            out = out,
            center = center,    # in boxlen unit
            width = width,      # in boxlen unit
            axis = axis,         # either ['x', 'y', 'z'] or a vector (list)
            field = field,
            weight_field = weight_field,
            kind = kind,       # either "slice" or "projection"
            normal = normal,        # None or vector (list)
            north = normal,         # None or vector (list)
            tag = tag,            # a free tag, used to distinguish
                                  # different versions where all other
                                  # parameters are the same
        )
        hashstr = dict_hash(params)
        h5fn = os.path.join(self.data_dir, hashstr + ".h5")
        if ds is None:
            ds = self.load_ds(out)
        if not os.path.exists(h5fn):
            if not os.path.isdir(self.data_dir):
                os.makedirs(self.data_dir)
            if kind == "slice":
                p = ds.slice(axis, field, center=center,
                            field_parameters={'width': width},
                            )
            elif kind == "projection":
                p = ds.proj(field=field, axis=axis, center=center,
                            weight_field=weight_field,
                            # weight_field=None,
                            field_parameters={'width': width},
                            )
            else:
                warning(f"Unknown kind: {kind}")
                return
            p.save_as_dataset(h5fn)
        data = yt.load(h5fn)
        if kind == "slice":
            p = yt.SlicePlot(data, axis, field)
        elif kind == "projection":
            p = yt.ProjectionPlot(data, axis, field,
                                  weight_field=weight_field)
        else:
            warning(f"Unknown kind: {kind}")
            return
        if field == ("gas", "density"):
            den_setup(p, weight_field=weight_field)
        if field == ("gas", "temperature"):
            T_setup(p)
        return p

    def projection_for_all_outputs(self, outdir, prefix="output",
                                   center='c', width=1.0):
        """Plot (density-weighted) projection of density for all frames of a
        simulation job.

        Args:
            ourdir (int):
            prefix (str): prefix of the filenames (default "output"). The name 
                of the figures would be prefix-x-output_00001.png
            center (str or list): (default 'c')
            width (float or tuple): (default 1.0) float as in boxlen unit or a
                tuple containing a number and unit, e.g. (0.1, 'pc')

        Returns:
            None

        """

        for i in self.get_all_outputs():
            for axis in ['x', 'y', 'z']:
                p = self.plot_prj(i, center=center, width=width, axis=axis)
                p.save(os.path.join(outdir, f"{prefix}-{axis}-output_{i:05d}.png"))
