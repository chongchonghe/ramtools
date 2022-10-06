#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot utilities for RAMSES outputs

Attributes
----------

Example
-------

"""

#-----------------------------------------------------------------------------
#    Author: Chong-Chong He
#    Data: 2021-02-09
#-----------------------------------------------------------------------------

import yt
from . import ytfast
from matplotlib.patches import Rectangle

def den_setup(p, zlim=None, time_offset=None, mu=1.4, unit='number_density',
              weight_field="density", is_time=False):
    """
    Args:
        p (YT plot):
        zlim (tuple): limits of the field
        mu (float): the mean particle weight. The mass density rho = mu * n
        quant (str): 'volume' or 'column' for volumetric and column density
    """

    if "density" in p.fields:
        fi = "density"
    elif ("gas", "density") in p.fields:
        fi = ("gas", "density")
    else:
        print("Failed to run den_setup: 'density' field not found")
        return
    if unit == 'number_density':
        if weight_field == "density" or weight_field == ("gas", "density"):
            p.set_unit(fi, 'cm**-3', equivalency='number_density',
                       equivalency_kwargs={'mu': mu})
    if zlim is not None:
        p.set_zlim(fi, zlim[0], zlim[1])
    # This will change the pixel size of the figure. Be careful!
    # p.set_figure_size(6)
    p.set_cmap(fi, 'inferno')
    if is_time:
        p.annotate_timestamp(time_format='{time:.3f} {units}',
                             time_offset=time_offset)

def T_setup(p, zlim=None, time_offset=None, is_time=False):
    if "temperature" in p.fields:
        fi = "temperature"
    elif ("gas", "temperature") in p.fields:
        fi = ("gas", "temperature")
    else:
        print("Failed to run T_setup: 'temperature' field not found")
        return
    if zlim is not None:
        p.set_zlim(fi, zlim[0], zlim[1])
    p.set_figure_size(6)
    p.set_cmap(fi, 'gist_heat')
    if is_time:
        p.annotate_timestamp(time_format='{time:.3f} {units}',
                             time_offset=time_offset)

# def T_setup(p, zlim=None, time_offset=None, is_time=True):
#     if zlim is not None:
#         p.set_zlim('temperature', zlim[0], zlim[1])
#     p.set_figure_size(6)
#     p.set_cmap('temperature', 'gist_heat')
#     if is_time:
#         p.annotate_timestamp(time_format='{time:.3f} {units}',
#                              time_offset=time_offset)

# def den_setup(p, zlim=None, time_offset=None, mu=1.4, unit='number_density'):
#     for field in ["density", ("gas", "density")]:
#         if unit == 'number_density':
#             p.set_unit(field, 'cm**-3',
#                        equivalency='number_density',
#                        equivalency_kwargs={'mu': mu})
#         if zlim is not None:
#             p.set_zlim(field, zlim[0], zlim[1])
#         # This will change the pixel size of the figure. Be careful!
#         p.set_figure_size(6)
#         p.set_cmap(field, 'inferno')
#         p.annotate_timestamp(time_format='{time:.3f} {units}',
#                              time_offset=time_offset)


def quick_plot_prj(ds, sinks, field='density', axis='z', use_h5=True,
                   **kwargs):
    """
    Args:
        ds: YT ds instance, e.g. ds = yt.load(...)
        sinks (n-by-3 array): position of the sink particles in code units
    """

    if use_h5:
        p = ytfast.ProjectionPlot(ds, axis, field,
                                  weight_field='density',
                                  **kwargs)
    else:
        p = yt.ProjectionPlot(ds, axis, field, weight_field=('gas', 'density'),
                              **kwargs)
    if field in ['density', ('gas', 'density')]:
        den_setup(p)
    # set the unit of the axes and of the sink positions to 'pc'
    p.set_axes_unit('pc')
    sinks = sinks / ds['boxlen']
    for pos in sinks:
        p.annotate_marker(pos, '.', coord_system='data',
                          plot_args={'color': 'cyan', 's': 40})
    return p


def plot_prj(ds, center, width, sinks, field='density', axis='x',
             axis_unit='pc', kind='prj', **kwargs):
    """
    Args:
        ds: YT ds instance, e.g. ds = yt.load(...)
        center: center in boxlen units (0, 1)
        width: width in boxlen units
        sinks (n-by-3 array): position of the sink particles in code units
    """

    if kind == 'prj':
        p = yt.ProjectionPlot(ds, axis, field, center=center, width=width,
                              weight_field='density', **kwargs)
    elif kind == 'slc':
        p = yt.SlicePlot(ds, axis, field, center=center, width=width, **kwargs)
    if field == 'density':
        den_setup(p)
    elif field == 'temperature':
        T_setup(p)
    # set the unit of the axes and of the sink positions to 'pc'
    p.set_axes_unit(axis_unit)
    sinks = sinks / ds['boxlen']
    for pos in sinks:
        # p.annotate_marker(pos, coord_system='data', plot_args={'color': 'cyan', 's': 40})
        p.annotate_marker(pos, coord_system='data', plot_args={'color': 'cyan'})
    return p


def overplot_time_tag(time, ax, loc='upper left', **kwargs):
    """
    Overplot time tag on top-left corner

    Args:
        time (float):
        ax:
        loc (str): one of 'upper left', 'upper right', 'upper center'
        timeshift:

    Returns:

    """
    if loc == 'upper left':
        pos = [.05, .95]
        va, ha = 'top', 'left'
    elif loc == 'upper right':
        pos = [.95, .95]
        va, ha = 'top', 'right'
    elif loc == 'upper center':
        pos = [.5, .95]
        va, ha = 'top', 'center'
    t = f"t = {time:.1f} Myr"
    ax.text(*pos, t, va=va, ha=ha, transform=ax.transAxes, **kwargs)


def add_scalebar(ax, length, label, h=0.014, left=0.03, right=None,
                 color='w', gap=0.01, **kwargs):
    """ Add a scalebar to a figure
    Author: ChongChong He

    Parameters
    ----------
    ax: matplotlib axes
        The axes passed to add_scalebar
    length: double
        The length of the scalebar in code units
    label: string
        The scalebar label
    h: double
        The height of the scalebar relative to figure height
    color: string
        color
    **kwargs: dict
        kwargs passed to ax.text

    Examples
    --------
    >>> im = plt.imread("/Users/chongchonghe/Pictures/bj.jpg")[:800, :800, :]
    >>> plt.imshow(im)
    >>> ax = plt.gca()
    >>> add_scalebar(ax, 200, '3 pc')
    >>> plt.show()

    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if right is None:
        left_pos = xlim[0] + left * (xlim[1] - xlim[0])
    else:
        left_pos = xlim[0] + right * (xlim[1] - xlim[0]) - length
    bottom = 0.03
    bottom_pos = ylim[0] + bottom * (ylim[1] - ylim[0])
    rect = Rectangle((left_pos, bottom_pos), length, h*(ylim[1] - ylim[0]),
                     facecolor=color, edgecolor=None)

    ax.add_patch(rect)
    ax.text(
        left_pos + 0.5*length,
        ylim[0] + (h + bottom + gap) * (ylim[1] - ylim[0]),
        label,
        ha='center', va='bottom',
        color=color,
        **kwargs,
        )


