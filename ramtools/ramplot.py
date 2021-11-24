# -*- coding: utf-8 -*-
"""
Some utilities to make figures.

Attributes:
    ...
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import yt
from typing import Dict, Any
import hashlib
import json

from ramtools import ytfast

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
                 use_h5=True,
                 force_redo=False,
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
            normal (tuple or list): normal vector
            north (tuple or lsit): north vector
            tag (any): a tag
            use_h5 (bool): whether or not to use h5 data
            force_redo (bool): whether or not to force remake the h5 data and figure

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
        if ds is None:
            ds = self.load_ds(out)

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

        if use_h5:
            if force_redo or (not os.path.exists(h5fn)):
                if not os.path.isdir(self.data_dir):
                    os.makedirs(self.data_dir)
                if kind == "slice":
                    if normal is None:
                        pass
                    p = ds.slice(axis, field, center=center,
                                field_parameters={'width': width},
                                )
                elif kind == "projection":
                    if normal is None:
                        p = ds.proj(field=field, axis=axis, center=center,
                                    weight_field=weight_field,
                                    # weight_field=None,
                                    field_parameters={'width': width},
                                    )
                    else:   # TODO
                        p = ds.proj(field=field, axis=axis, center=center,
                                    weight_field=weight_field,
                                    # weight_field=None,
                                    field_parameters={'width': width},
                                    )
                else:
                    warning(f"Unknown kind: {kind}")
                    return
                p.save_as_dataset(h5fn, ) #fields=[field])
                # # option 2, FixedResolutionBuffer
                # if axis in [2, 'z']:
                #     bounds = (center[0] - width/2, center[1] - width/2,
                #               center[1] - width/2, center[1] + width/2)
                # frb = yt.FixedResolutionBuffer(p, bounds, (1024, 1024))
                # frb.save_as_dataset(h5fn, fields=[field])
                # # frb.export_hdf5(h5fn, fields=[field])

            # data = yt.load(h5fn).all_data()
            data = yt.load(h5fn)
            if kind == "slice":
                p = yt.SlicePlot(data, axis, field, center=center)
            elif kind == "projection":
                p = yt.ProjectionPlot(data, axis, field, center=center,
                                      width=width, weight_field=weight_field)
            else:
                warning(f"Unknown kind: {kind}")
                return
        else:       # not using h5 data
            if kind == "slice":
                p = yt.SlicePlot(data, axis, field)
            elif kind == "projection":
                p = yt.ProjectionPlot(ds, axis, field, center=center,
                                      width=width, weight_field=weight_field)
        if field == ("gas", "density"):
            den_setup(p, weight_field=weight_field)
        if field == ("gas", "temperature"):
            T_setup(p)
        # del data # cleanup memory
        return p

    def projection_for_all_outputs(self, outdir, prefix="output",
                                   center='c', width=1.0, force_redo=False):
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

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for i in self.get_all_outputs():
        # for i in range(19, self.get_all_outputs()[-1]+1):
            for axis in ['x', 'y', 'z']:
                fn = os.path.join(outdir, f"{prefix}-{axis}-output_{i:05d}.png")
                print("\nPlotting", fn)
                try:
                    p = ytfast.ProfilePlot(
                        ds = self.load_ds(i),
                        axis = axis,
                        fields = ('gas', 'density'),
                        center = 'c',
                        width = 0.8,
                        axes_unit = 'pc',
                        weight_field = ('gas', 'density'),
                    )
                    # p = self.plot_prj(i, center=center, width=width, axis=axis,
                    #                   force_redo=force_redo)
                except Exception as _e:
                    print("Caught error (ignored):")
                    print(_e)
                    continue
                p.save(fn)
                print("Done")
                del p

    def plot_2d_profile(self, out, center, radius,
                        field_x, field_y, field, ds=None, weight_field=None,
                        lims={}, units={}, vlims=None, with_cb=True,
                        force_redo=False, define_field=None):
        """ Plot 2D profile. Will store intermediate to disk for fast recreate of
        the figure.

        Args:
            jobpath (str): path to job
            out (int):
            center (tuple or str): tuple (boxlen unit) or 'c' (center)
            radius (float or tuple): radius (half-width) of the region
            ds (YT ds): optional, feed ds (default None)
            lims (dict): limits of the axes, e.g. {'density': [1e3, 1e5],
                'temperature': [1e0, 1e4]}
            vlims (tuple): vlims of the colorbar (default None)
            with_cb (bool): toggle colorbar (default True)
            define_field (function): a function that defines a new field. It
                should take ds as an argument.

        Returns:
            yt profile plot

        """

        from matplotlib import colors

        if not os.path.exists(self.get_info_path(out)):
            return None, None
        h5dir = os.path.join(self.data_dir, "profiles")
        if not os.path.exists(h5dir):
            os.makedirs(h5dir)
        radius = self.to_boxlen(radius)
        left = tuple([c - radius for c in center])
        right = tuple([c + radius for c in center])
        if isinstance(center, list):
            center = tuple(center)
        n_bins=(128, 128)
        logs = {field_x: True, field_y: True}
        lims_str = {str(key): value for key, value in lims.items()}
        logs_str = {str(key): value for key, value in logs.items()}
        units_str = {str(key): value for key, value in units.items()}
        hash_params = dict(
            jobdir = os.path.abspath(self.jobPath),
            out = out, center = center, radius = radius,
            lims = lims_str, n_bins = n_bins, logs = logs_str,
            units = units_str
        )
        # print(hash_params)
        hashstr = dict_hash(hash_params)
        h5fn = os.path.join(h5dir, hashstr + ".h5")
        if force_redo or (not os.path.exists(h5fn)):
            if ds is None:
                ds = self.load_ds(out)
            if define_field is not None:
                define_field(ds)
            box = ds.box(left, right)
            p = yt.create_profile(box, [field_x, field_y], field,
                                  weight_field=weight_field, logs=logs,
                                  extrema=lims, n_bins=n_bins,
                                  units=units,
                                 )
            p.save_as_dataset(h5fn)
        prof = yt.load(h5fn)
        dat = prof.data[field].T
        extents = [prof.data[field_x].min(), prof.data[field_x].max(),
                   prof.data[field_y].min(), prof.data[field_y].max()]
        extent_log = np.log10(extents)
        if vlims is None:
            themax = dat.max()
            if not dat.min() > 0.0:
                themin = themax * 1e-6
            else:
                themin = dat.min()
            vlims = [themin, themax]
        f, ax = plt.subplots()
        ax.imshow(dat,
                  norm=colors.LogNorm(vmin=vlims[0], vmax=vlims[1]),
                  extent=extent_log,
                  aspect="auto", origin="lower",
                  )
        str1 = field_x[1] if isinstance(field_x, tuple) else field_x
        str2 = field_y[1] if isinstance(field_y, tuple) else field_y
        ax.set(xlabel="log " + str1, ylabel=r"log " + str2)
        del dat
        try:
            del p
        except:
            pass
        return f, ax

    def plot_phaseplot_for_all(self, figdir, center, radius, prefix='', phaseplot_kwargs={}):
        """

        Usage:

            >>> plot_phaseplot_for_all(center='c', radius=0.4, x_field=('gas', 'density'),
                    y_field=('gas', 'temperature'), z_fields=('gas', 'cell_mass'))
        """

        if not os.path.exists(figdir):
            os.makedirs(figdir)
        for out in self.get_all_outputs():
            ds = self.load_ds(out)
            sph = yt.Sphere(center=center, radius=radius, )
            f, ax, p, cb = ytfast.PhasePlot(sph, **phaseplot_kwargs)
            figfn = f"{figdir}/{prefix}_{out:03d}.pdf"
            f.savefig(figfn)
