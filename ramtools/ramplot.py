# -*- coding: utf-8 -*-
"""
Some utilities to make figures.

Attributes:
    ...
"""

import os
import logging
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import yt
from typing import Dict, Any
import hashlib
import json

# from ramtools import ytfast
from . import ytfast, utilities
from .units import AU
from .plotutils import T_setup, den_setup

try:
    yt.set_log_level(40)
except:
    pass

logger = logging.getLogger()
logger.setLevel(logging.INFO)
# logger.debug("Debugging is turned on")
# logger.error("Error is turned on")

from .ramses import Ramses

def dict_hash(dictionary: Dict[str, Any]) -> str:
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # Sort arguments so {'a': 1, 'b': 2} is the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()



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
                    # p = ds.slice(axis, field, center=center,
                    #             field_parameters={'width': width},
                    #             )
                    p = yt.SlicePlot(ds, axis, field, center=center,
                                     width=width).data_source
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
                print(h5fn, 'saved.')
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

        ytfast.set_data_dir(os.path.join(self.data_dir, 'profiles'))
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        for out in self.get_all_outputs():
            if not out % 20 == 0:
                continue
            ds = self.load_ds(out)
            sph = ds.sphere(center=center, radius=radius, )
            timeit.a()
            f, ax, p, cb = ytfast.PhasePlot(sph, **phaseplot_kwargs)
            timeit.b()
            figfn = f"{figdir}/{prefix}_{out:03d}.png"
            f.savefig(figfn, dpi=300)



def to_boxlen(quant, ds):
    if type(quant) == tuple:
        return float(yt.YTQuantity(*quant) / ds.length_unit)
    else:
        return quant


def plot_a_region(
        ram, out, ds, center,
        fields='den', kind='slc', axis='z', width=None,
        center_vel=None, L=None, direcs='face',
        zlims={}, l_max=None, is_id=False,
        bar_length=2e3, sketch=False,
        time_offset='tRelax',
        is_set_size=True,
        is_time=True,
        streamplot_kwargs={},
        more_kwargs={},
):
    """Do projection or slice plots in a region.

    Args:
        ram: a ramtools.Ramses instance
        out (int): output frame
        ds (yt.ds): yt.load instance
        center (list_like): the center of the region like in yt.SlicePlot
        fields (str or tuple): the field to plot. The avaialble options are: 'den' or 'density' - density. 'logden' - log of density. 'T' or 'temperature' - temperature. 'T_vel_rela' - temperature overplot with velocity field. 'pressure' - thermal pressure. 'magstream' - stream lines of magnetic fields on top of density slice. The magnetic strength is be indicated by the color of the stream lines. 'beta' - plasma beta parameter. 'AlfMach' - Alfvenic Mach number. 'xHII' - hydrogen ionization fraction. 'mach' - thermal Mach number. 'mag' - magnetic field strength. 'vel' - velocity field on top of density slice. 'vel_rela' - relative velocity field, i.e. the velocity field in the frame with a velocity defined by center_vel. 'mach2d' - thermal Mach number in the plane. 'p_mag' - magnetic pressure. 
        kind (str): 'slc' or 'prj'. Default: 'slc'
        axis (str or int): One of (0, 1, 2, 'x', 'y', 'z').
        width (float or tuple): width of the field of view. e.g. 0.01,
            (1000, 'AU').
        center_vel (list_like or None): the velocity of the center. (Default
            None)
        L (list_like): the line-of-sight vector. Will overwrite axis.
            Default: None
        direcs (str): 'face' or 'edge'. Will only be used if L is not None.
        zlims (dict): The limits of the fields. e.g. {'density': [1e4, 1e8]},
            {'B': [1e-5, 1e-2]} (the default of B field).
        l_max (int):
        is_id (bool): toggle marking sink ID (the order of creation)
        bar_length (float): not using
        sketch (bool): If True, will do a faster plot of a simple imshow for
            test purpose. Default: False.
        is_time (bool): toggle overplotting time tag
        time_offset (str or int or float): One of the following cases:
            str: 'rRelax'
            int: the frame of which the time is used
            float: time in Myr
        is_set_size (bool): toggle setting figure size to 6 to make texts
            bigger. Default: True.
        streamplot_kwargs (dict): More kwargs for plt.streamplot that is used when fields='magstream'. Default: {}
        more_kwargs (dict): more field-specific kwargs. Default: {}

    Returns:
        yt figure or plt.figure

    """

    from matplotlib import colors

    r = ram
    ytfast.set_data_dir(os.path.join(r.ramses_dir, 'h5_data'))
    assert width is not None
    width = to_boxlen(width, r.ds1)
    if isinstance(time_offset, int):
        time_offset = r.get_time(time_offset)
    elif isinstance(time_offset, str):
        if time_offset is 'tRelax':
            time_offset = r.tRelax
    else:
        time_offset = time_offset
    if l_max is None: l_max = ds.max_level
    use_h5 = False      # TODO: enable use_h5

    # calculated north and right, if L is not None
    if L is None:
        if axis in [0, 1, 2]:
            axis = ['x', 'y', 'z'][axis]
        is_onaxis = True
        los = axis
        axismap = {'x': [1, 0, 0], 'y': [0, 1, 0], 'z': [0, 0, 1]}
        los_v = axismap[axis]
        north = axismap[{'x': 'z', 'y': 'x', 'z': 'y'}[axis]]
        logger.debug('is_onaxis = {}'.format(is_onaxis))
    else:
        is_onaxis = False
        L /= norm(L)
        tmpright = np.cross([0, 0, 1], L)
        tmpright = tmpright / norm(tmpright)
        if direcs == 'face':
            los = L
            north = -1 * tmpright
        elif direcs == 'edge':
            los = tmpright
            north = L
        los_v = los
    right = np.cross(north, los_v)
    logger.debug('is_onaxis = {}'.format(is_onaxis))
    logger.info(f'doing {r.jobPath}, los={los}, north={north}, right={right}')

    # ------------------ doing the plot --------------------
    if 'density' in zlims:
        den_zlim = zlims['density']
    elif 'den' in zlims:
        den_zlim = zlims['den']
    else:
        den_zlim = [None, None]
    gas_field = fields
    is_log = True
    zlim = None
    cb_label = None
    cmap = 'viridis'
    if fields in ['T', 'temperature', 'T_vel_rela']:
        gas_field = 'temperature'
    if fields in ['pressure', 'p_mag']:
        cmap = 'gist_heat'
    if fields in ['beta', 'magbeta', 'AlfMach']:
        cmap = 'seismic'
        add_beta_corrected(ds)
        zlim = [-3, 3]
        is_log = False
        if fields == 'magbeta':
            gas_field = "logPlasmaBeta"
            def _logBeta(field, data):
                return np.log10(data['beta'])
            ds.add_field(("gas", gas_field), function=_logBeta)
            cb_label = "log (Plasma Beta)"
        if fields == 'AlfMach':
            gas_field = ("gas", "LogMachA")
            def _LogMachA(field, data):
                """ MachA = v_rms / v_A = gamma / 2 * Mach^2 * beta """
                gamma = 5 / 3
                return np.log10(gamma / 2 * data["mach"] ** 2 * data["beta"])
            ds.add_field(gas_field, function=_LogMachA)
            cb_label = "log $M_A$"
    if fields in ['xHII']:
        zlim = [1e-6, 1]
    if fields in ['mach', 'AlfMach']:
        add_mach(ds, center_vel)
        if gas_field == 'mach':
            gas_field = 'logmach'
            add_logmach(ds)
            cb_label = r'log ($\mathcal{M}_s$)'
            is_log = False
            cmap = 'seismic'
            zlim = [-2, 2]
    if fields in ['den', 'density', 'mag', 'vel', 'vel_rela']:
        gas_field = 'density'
    if fields in ['logden', ('gas', 'logden')]:
        def _logden(_field, data):
            mu = 1.4
            return data[('gas', 'density')] / (mu * yt.physical_constants.mass_hydrogen)
        gas_field = ('gas', 'logden')
        ds.add_field(gas_field, function=_logden, take_log=True,
                     sampling_type='cell', units="cm**(-3)")
    if fields == 'magstream':
        def _rel_B_1(_field, data):
            Bx = (data["x-Bfield-left"] + data["x-Bfield-right"]) / 2.
            By = (data["y-Bfield-left"] + data["y-Bfield-right"]) / 2.
            Bz = (data["z-Bfield-left"] + data["z-Bfield-right"]) / 2.
            B = yt.YTArray([Bx, By, Bz])
            return np.tensordot(right, B, 1)
        def _rel_B_2(_field, data):
            Bx = (data["x-Bfield-left"] + data["x-Bfield-right"]) / 2.
            By = (data["y-Bfield-left"] + data["y-Bfield-right"]) / 2.
            Bz = (data["z-Bfield-left"] + data["z-Bfield-right"]) / 2.
            B = yt.YTArray([Bx, By, Bz])
            return np.tensordot(north, B, 1)
        ds.add_field(('gas', 'rel_B_1'), function=_rel_B_1)
        ds.add_field(('gas', 'rel_B_2'), function=_rel_B_2)
    if fields in ['vel_rela', 'mach2d']:
        assert north is not None
        def _rel_vel_1(_field, data):
            rel_vel = yt.YTArray(
                [data['velocity_x'] - yt.YTArray(center_vel[0], 'cm/s'),
                 data['velocity_y'] - yt.YTArray(center_vel[1], 'cm/s'),
                 data['velocity_z'] - yt.YTArray(center_vel[2], 'cm/s')]
            )
            return yt.YTArray(np.tensordot(right, rel_vel, 1), 'cm/s')
        def _rel_vel_2(_field, data):
            rel_vel = yt.YTArray(
                [data['velocity_x'] - yt.YTArray(center_vel[0], 'cm/s'),
                 data['velocity_y'] - yt.YTArray(center_vel[1], 'cm/s'),
                 data['velocity_z'] - yt.YTArray(center_vel[2], 'cm/s')]
            )
            return yt.YTArray(np.tensordot(north, rel_vel, 1), 'cm/s')
        def _mach2d(_field, data):
            rel_speed_sq = data['rel_vel_1'].in_cgs().value ** 2 + \
                           data['rel_vel_2'].in_cgs().value ** 2  # cm/s
            T = data['temperature'].value
            # km/s, assuming gamma=5/3. mu is inside T, therefore this
            # is precise.
            cs = 0.11729 * np.sqrt(T)
            return 1e-5 * np.sqrt(rel_speed_sq) / cs  # dimensionless
        ds.add_field(('gas', 'rel_vel_1'),
                     function=_rel_vel_1, units="cm/s")
        ds.add_field(('gas', 'rel_vel_2'),
                     function=_rel_vel_2, units="cm/s")
        ds.add_field(('gas', 'mach2d'), function=_mach2d, )
        if fields == 'mach2d':
            gas_field = 'logmach2d'
            def _logmach2d(_field, data):
                return np.log10(data["mach2d"])
            ds.add_field(('gas', gas_field), function=_logmach2d, )
            cb_label = r'log ($\mathcal{M}_{s, 2D}$)'
            is_log = False
            cmap = "seismic"
            zlim = [-2, 2]

    if fields != 'magstream':
        if kind == 'slc':
            if is_onaxis:
                p = yt.SlicePlot(ds, los, gas_field, width=width,
                                 center=center)
            else:
                p = yt.OffAxisSlicePlot(
                    ds, los, gas_field, width=width,
                    center=center, north_vector=north)
        elif kind == 'prj':
            if is_onaxis:
                p = ytfast.ProjectionPlot(
                    ds, los, gas_field, width=width,
                    center=center, weight_field='density')
            else:
                print("Doing OffAxis projection plot...")
                p = yt.OffAxisProjectionPlot(
                    ds, los, gas_field, width=width,
                    center=center, weight_field='density',
                    max_level=l_max,
                    north_vector=north)
                print("Done")
        elif kind == 'colden':
            if is_onaxis:
                p = ytfast.ProjectionPlot(
                    ds, los, gas_field, width=width,
                    center=center, weight_field=None)
            else:
                p = yt.OffAxisProjectionPlot(
                    ds, los, gas_field, width=width,
                    center=center, weight_field=None,
                    max_level=l_max,
                    north_vector=north)
        # p.set_colorbar_label(field, _label)
        if is_set_size:
            p.set_figure_size(6)
        p.set_axes_unit('AU')
        if gas_field in zlims:
            if zlims[gas_field] is not None:
                p.set_zlim(gas_field, *zlims[gas_field])
        if is_time:
            p.annotate_timestamp(time_format='{time:.2f} {units}',
                                 time_offset=time_offset,
                                 text_args={'fontsize':8,
                                            'color':'w'},
                                 corner='upper_left',
                                 )
        # Overplot sink particles
        args = {}
        if 'sink_colors' in more_kwargs.keys():
            args['colors'] = more_kwargs['sink_colors']
        if 'mass_lims' in more_kwargs.keys():
            args['lims'] = more_kwargs['mass_lims']
        if 'colors' in more_kwargs.keys():
            args['colors'] = more_kwargs['colors']
        r.overplot_sink_with_id(
            p, out, center, width/2, is_id=is_id,
            zorder='mass', withedge=1, **args)

        # set cmap, zlim, and other styles
        p.set_log(gas_field, is_log)
        p.set_cmap(gas_field, cmap)
        if zlim is not None:
            p.set_zlim(gas_field, zlim)
        p.set_colorbar_label(gas_field, cb_label)
        if gas_field == 'density':
            den_setup(p, den_zlim, time_offset=time_offset)
        if gas_field == 'temperature':
            T_setup(p, [3, 1e4], time_offset=time_offset)

        if fields == 'vel':
            p.annotate_velocity(factor=30, scale=scale,
                                scale_units='x',
                                plot_args={"color": "cyan"})
            p.set_colorbar_label(gas_field, r'log ($\mathcal{M}_s$)')
        if fields == ['vel_rela', 'T_vel_rela']:
            p.annotate_quiver('rel_vel_1', 'rel_vel_2', factor=30, scale=scale,
                              scale_units='x',
                              plot_args={"color": "cyan"})
        if fields in ['vel', 'vel_rela', 'T_vel_rela']:
            # scale ruler for velocities
            ruler_v_kms = 4     # km/s
            width_au = width * r.boxlen * 2e5
            bar_length = 2e3 * width_au / 2e4
            coeff = bar_length            # length of ruler in axis units [AU]
            scale = 1e5 * ruler_v_kms / coeff  # number of cm/s per axis units
            # add a scale ruler
            p.annotate_scale(coeff=coeff, unit='AU',
                             scale_text_format="{} km/s".format(ruler_v_kms))
        return p
    else:           # streamplot
        if not sketch:
            sl = ds.cutting(los_v, center, north_vector=north)
            hw = width / 2.
            bounds = [-hw, hw, -hw, hw]
            size = 2 ** 10
            frb = yt.FixedResolutionBuffer(sl, bounds, [size, size])
            m_h = 1.6735575e-24
            mu = 1.4
            den = frb['density'].value / (mu * m_h)
            unitB = utilities.get_unit_B(ds)
            logger.info(f"unitB.value = {unitB.value}")
            Bx = frb['rel_B_1'].value * unitB.value  # Gauss
            By = frb['rel_B_2'].value * unitB.value  # Gauss
            Bmag = np.sqrt(Bx ** 2 + By ** 2)
            logger.info(f"Bmag min = {Bmag.min()}, max = {Bmag.max()}")
            hw_au = hw * r.unit_l / AU
            grid_base = np.linspace(-hw_au, hw_au, size)
            bounds_au = [-hw_au, hw_au, -hw_au, hw_au]
        else:
            den = np.random.random([3, 3])
            hw_au = 10
            den_zlim = [None, None]
            bounds_au = [-hw_au, hw_au, -hw_au, hw_au]

        # fig, ax = plt.subplots()
        if 'figax' in more_kwargs.keys():
            fig, ax = more_kwargs['figax']
        else:
            fig, ax = plt.subplots()
        im = ax.imshow(np.log10(den), cmap='inferno', extent=bounds_au,
                       vmin=den_zlim[0], vmax=den_zlim[1], origin='lower')
        colornorm_den = colors.Normalize(vmin=den_zlim[0], vmax=den_zlim[1])
        ax.set(
            xlabel="Image x (AU)",
            ylabel="Image y (AU)",
            xlim=[-hw_au, hw_au],
            ylim=[-hw_au, hw_au]
        )
        cmap = mpl.cm.get_cmap("Greens", 6, )
        Bvmin, Bvmax = -5, -2
        if "Blims_log" in more_kwargs.keys():
            Bvmin, Bvmax = more_kwargs["Blims_log"]
        if "B" in zlims.keys():
            Bvmin = np.log10(zlims['B'][0])
            Bvmax = np.log10(zlims['B'][1])
        colornorm_stream = colors.Normalize(vmin=Bvmin, vmax=Bvmax)
        # streamplot_kwargs = {}
        # if "streamplot" in more_kwargs.keys():
        #     streamplot_kwargs = more_kwargs["streamplot"]
        # linewidth = 0.2
        if "stream_linewidth" in more_kwargs.keys():
            linewidth = more_kwargs["stream_linewidth"]
        if not sketch:
            logmag = np.log10(Bmag)
            # scaled_mag = (logmag - Bvmin) / (Bvmax - Bvmin)
            strm = ax.streamplot(grid_base, grid_base, Bx, By,
                                 color=np.log10(Bmag),
                                 density=1.2,
                                 # linewidth=linewidth,
                                 # linewidth=1.5 * scaled_mag,
                                 # cmap='Greens',
                                 cmap=cmap,
                                 norm=colornorm_stream,
                                 arrowsize=0.5,
                                 **streamplot_kwargs,
                                 )
        # plt.subplots_adjust(right=0.8)
        plot_cb = True
        plot_cb2 = True
        if 'plot_cb' in more_kwargs.keys():
            plot_cb = more_kwargs['plot_cb']
        if 'plot_cb2' in more_kwargs.keys():
            plot_cb2 = more_kwargs['plot_cb2']
        if plot_cb:
            if 'cb_axis' in more_kwargs.keys():
                ax2 = more_kwargs['cb_axis']
                # pos0 = ax.get_position()
                # cbaxis = fig.add_axes([pos0.x1, pos0.y0, 0.06, pos0.height])
                # cb = plt.colorbar(im, cax=cbaxis)
                cb2 = mpl.colorbar.ColorbarBase(
                    ax2, orientation='vertical',
                    cmap=mpl.cm.get_cmap("inferno"),
                    norm=colornorm_den,
                )
            else:
                cb2 = fig.colorbar(im, ax=ax, pad=0)
            cb2.set_label("log $n$ (cm$^{-3}$)")
        if plot_cb2:
            if 'cb2_axis' in more_kwargs.keys():
                ax2 = more_kwargs['cb2_axis']
            elif plot_cb:
                pos1 = ax.get_position()
                ax2 = fig.add_axes([pos1.x1 + .15, pos1.y0, 0.02, pos1.height])
            else:
                pos1 = ax.get_position()
                ax2 = fig.add_axes([pos1.x0, pos1.y0, 0.02, pos1.height])
            colornorm2 = colors.Normalize(vmin=den_zlim[0], vmax=den_zlim[1])
            cb2 = mpl.colorbar.ColorbarBase(
                ax2, orientation='vertical', cmap=cmap, norm=colornorm_stream,
                # ticks=[-4, -3, -2]
            )
            cb2.set_label("log B (Gauss)")
        # plt.savefig(fs['magstream'], dpi=300)
        return fig

