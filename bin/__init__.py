#!/usr/bin/env python

import os, sys
from math import log10
import argparse
import warnings
warnings.filterwarnings("ignore")

# put imports inside main() to make -h run faster
import yt
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from ramtools import ytfast, Ramses, plotutils, ramses, utilities as ut
from matplotlib.pyplot import cm

try:
    import scienceplots
    plt.style.use(['science', 'no-latex'])
except:
    pass


def parse_outs(outs_str):
    outs = []
    for par in outs_str.split(','):
        if '-' in par:
            splits = [int(i) for i in par.split('-')]
            if len(splits) == 2:
                outi, outf = splits
                outs += list(range(outi, outf + 1))
            elif len(splits) == 3:
                outi, outf, diff = splits
                outs += list(range(outi, outf + 1, diff))
        else:
            outs += [int(par)]
    return outs


def get_args_prj():
    args = argparse.ArgumentParser(description=("""
A program to efficiently make slices/projections plots of a RAMSES simualtion. 

Quick examples:
>>> plot_projection ..
>>> plot_projection . 10-20 z -t figures
>>> plot_projection . 10-20 x,y -t figures -c 0.6,0.6,0.6 -w 10,pc \
-f den,temp -k slc --sink --den-zlim 1e4 1e10 --T-zlim 10 1e4
"""),
        formatter_class=argparse.RawTextHelpFormatter)
    args.add_argument("jobdir", nargs="?", default="..",
                      help="(Default '..') Path to the simulation directory")
    args.add_argument("outs",  nargs="?", default="all",
                      help="(Default 'all') Output frames. If 'all', do all "
                           "output frames. "
                      "Examples: '10', '10,11,13', '10-20', '10-20-2'")
    args.add_argument("axes", nargs="?", default="x,y,z",
                      help="(Default 'x,y,z') The line of sight. Examples: "
                           "'x', 'x,y,z'")
    args.add_argument("-c", "--center", default='c',
                      help="(Default 'c') The center.\n"
                           "Cases:\n"
                           "    'c': = [0.5, 0.5, 0.5]\n"
                           "    'peak': the peak density location"
                           "    x,y,z: position in boxlen unit, e.g. 0.5,0.5,0.5\n"
                           "    a single integer: use the location of a sink "
                           "    particle with the given index as the center "
                           "    (starting from index 0). You may use --inzoombox "
                           "    to restrict sinks in zoom box only.\n"
                     )
    args.add_argument("-l", "--los", type=float, nargs=3,
                      help="The line of sight vector, usually the face-on "
                           "vector of the disk. When this is set, will plot a "
                           "face-on and edge view that are perpendicular to "
                           "each other."
                      )
    args.add_argument("--no-edgeon", action="store_true",
                      help="Turn off plotting edge-on view when los is set")
    args.add_argument("-w", "--width", default='1',
                      help="The width of the view. Examples: '1', '1000,AU'. Default: '1'")
    args.add_argument("-f", "--fields", default='den',
                      help="The fields to plot. Examples: 'den', 'den,temp'. Default: 'den'")
    args.add_argument("-t", "--to", default='./projections',
                      help="The destination directory. Default: "
                           "'./projections'")
    args.add_argument("-k", "--kinds", default='slc,prj',
                      help="Plot types, either 'slc', 'prj', or 'slc,"
                           "prj' (default)")
    args.add_argument("-s", "--sink", action="store_true",
                      help="Turn on overplotting sink particles")
    args.add_argument("--sink-color", help="Set sink particle color. Default: Greens. Examples: Greens_r, Reds, Reds_r")
    args.add_argument("--mass-lim", type=float, nargs=2,
                      help="Colormap limits (in Msun) of the sink particles")
    args.add_argument("-o", "--overwrite", action="store_true",
                      help="Turn on overwritting existing figures. If left "
                           "off (default), will skip existing files.")
    args.add_argument("--overwrite-cache", action="store_true",
                      help="Turn on overwritting cache files used by ytfast")
    args.add_argument("-a", "--annotate-scale", action="store_true",
                      help="Toggle annotate scale")
    args.add_argument("--show-axes", action="store_true",
                      help="Toggle show axes")
    args.add_argument("--scale", help="Coeff and unit of annotate_scale")
    args.add_argument("--scaleloc", help="Location of the scale bar: lower_right (default) or"
                      "lower_left")
    args.add_argument("--no-ytfast", action="store_true",
                      help="Turn off using ytfast to speedup projection/slice plot")
    args.add_argument("--den-zlim", type=float, nargs=2,
                      help="zlim of density (in cm-3)")
    args.add_argument("--T-zlim", type=float, nargs=2,
                      help="zlim of temperature (in K)")
    args.add_argument("--xHII-zlim", type=float, nargs=2,
                      help="zlim of xHII")
    args.add_argument("--contour", type=float, help="column density threshold for contour")
    args.add_argument("--cmap", type=str, help="colormap")
    args.add_argument("--dry-run", action="store_true",
                      help="Toggle dry run: making one figure only and store it in local dierctory")
    args.add_argument("--inzoombox", action="store_true",
                      help="with center=integer or center='peak', use sink particles/gas in zoombox only.")
    args.add_argument("--rise", help="raise the slice plot by percent of the box width")
    args.add_argument("--t0", help="set the t0 and annotate time. t0 could be one of "
                                   "the following options: 1. an integer indicating the "
                                   "time of the corresponding output frame; 2. a float "
                                   "indicating the t0 in Myr")
    args.add_argument("--time-unit", help="time unit, either Myr (default) or kyr")
    args.add_argument("--darkbg", action="store_true", help="Use dark background")
    args.add_argument("--draw-box", help="If not None, draw a box with this width")
    args.add_argument("--box-center", help="define the center of the box to draw, used along with --draw-box. "
                                           "The syntax is similar to center")
    args.add_argument("--skip", help="skip certain outputs")
    # LEGACY parameters
    args.add_argument("--hide-axes", action="store_true",
                      help="LEGACY. Now as the default. Use --show-axes to show axes")
    args.add_argument("--prefix", default="", 
                      help="add prefix to image filename")
    args.add_argument("--figsize", type=float, default=6.0, help="figure size")
    args.add_argument("--den-limit-lower", type=float, help="cut off the density below the limit")
    return args.parse_args()


def plot_projection(args):
    """
    Main function of plot_projection

    """

    # # put imports inside main() to make -h run faster
    # import yt
    # import numpy as np
    # import matplotlib.pyplot as plt
    # from glob import glob
    # from ramtools import ytfast, Ramses, plotutils, ramses, utilities as ut
    # from matplotlib.pyplot import cm
    # import scienceplots

    # plt.style.use(['science', 'no-latex'])

    if args.darkbg:
        plt.style.use('dark_background')

    r = Ramses(args.jobdir)
    zoomc = r.read_zoom_center()
    zoomr = r.read_zoom_radius()
    if args.outs == 'all':
        outs = r.get_all_outputs()
    else:
        outs = parse_outs(args.outs)
        # outs = []
        # for par in args.outs.split(','):
        #     if '-' in par:
        #         splits = [int(i) for i in par.split('-')]
        #         if len(splits) == 2:
        #             outi, outf = splits
        #             outs += list(range(outi, outf + 1))
        #         elif len(splits) == 3:
        #             outi, outf, diff = splits
        #             outs += list(range(outi, outf + 1, diff))
        #     else:
        #         outs += [int(par)]
    if args.rise is not None:
        rises = parse_outs(args.rise)
    # use_sink_as_center = False
    # use_gas_peak_as_center = False

    if 'zoomr' in args.width:
        nmls = sorted(glob(f"{args.jobdir}/*.nml"))
        if len(nmls) == 0:
            print(f"Failed to load a namelist file. Is {args.jobdir} a valid job directory?")
            return
        nml = nmls[0]
        # print(f"Reading width from namelist file {nml}")
        idx = -1
        if len(args.width) > 5:
            idx = int(args.width[5:])
        width = ut.read_zoom_radius(nml, idx) * 2
    elif ',' in args.width:
        p1, p2 = args.width.split(',')
        width = (float(p1), p2)
    else:
        width = float(args.width)
    if args.draw_box is not None:
        if ',' in args.draw_box:
            p1, p2 = args.draw_box.split(',')
            box_width = (float(p1), p2)
        else:
            box_width = float(args.draw_box)
    os.makedirs(args.to, exist_ok=1)
    # time offset
    if args.t0 is not None:
        try:
            timeoffset = (r.get_time(int(args.t0)), 'Myr')
        except ValueError:
            timeoffset = (float(args.t0), 'Myr')
    axes = None
    if args.los is None:
        axes = args.axes.split(',')
    else:
        if args.no_edgeon:
            axes = ["face"]
        else:
            axes = ["face", "edge"]
        face = np.array(args.los)
        right = np.cross([0, 0, 1], face)
        left = -1. * right
        dires = {"face": [face, left],
                 "edge": [right, face]}
    skips = []
    if args.skip is not None:
        for i in args.skip.split(','):
            skips.append(int(i))
    for out in outs:
        if not r.exist(out) or out in skips:
            print(f"Skipping frame {out}")
            continue
        ds = r.load_ds(out)
        width = ramses.to_boxlen(width, ds)
        if args.draw_box is not None:
            box_width = ramses.to_boxlen(box_width, ds)

        # if args.den_limit_lower is not None, cut off the density below the limit
        if args.den_limit_lower is not None:
            ad = ds.all_data()
            high_density_region = ad.cut_region([f"(obj[('gas', 'density')] > {args.den_limit_lower})"])

        def parse_center(_center):
            if _center == 'c':
                return [0.5, 0.5, 0.5]
            elif _center == 'zoomc':    # use zoom center
                return zoomc
            elif _center == 'peak':    # use peak density location
                # use_gas_peak_as_center
                print(f"Finding the peak density location in out {out}...")
                # zoomr3 = r.read_zoom_radius(-2)
                zoomr3 = zoomr
                return plotutils.find_peak_density_location(
                    ds,
                    [izoomc - zoomr3 for izoomc in zoomc],
                    [izoomc + zoomr3 for izoomc in zoomc],
                )
            elif ',' in _center:
                return [float(i) for i in _center.split(',')]
            else:   # _center is int, use sink particle location
                sink_center = int(_center)
                try:
                    pos = r.get_sink_positions(out) / r.boxlen
                    if args.inzoombox:
                        inside = np.max(pos - zoomc, axis=1) <= zoomr
                        if np.sum(inside) == 0:
                            print(f"No sink particles inside zoombox in frame {out}")
                            return
                        return pos[inside][sink_center]
                    else:
                        return pos[sink_center]
                    # center = r.get_sink_positions(out)[sink_center] / r.boxlen
                except FileNotFoundError:
                    return
                except ramses.NoSinkParticle:
                    print(f"No sink particle found in out {out}.")
                    return

        center = parse_center(args.center)
        if center is None:
            print(f"Could not find a proper center in out {out}.")
            continue
        if args.draw_box is not None:
            boxcenter = parse_center(args.box_center)
            if boxcenter is None:
                print(f"Could not find a proper boxcenter in out {out}.")
                continue

        for field in args.fields.split(','):
            if field == "den":
                thefield = ("gas", "density")
            elif field in ["tem", "temp"]:
                thefield = ("gas", "temperature")
            elif field == "grid_level":
                thefield = ("index", field)
            else:
                thefield = ("gas", field)
            for kind in args.kinds.split(','):
                for axis in axes:

                    # rise slices
                    if args.rise and kind == 'slc':
                        assert axis in 'xyz'
                        if axis == 'x':
                            centers = [[center[0] + width * (i - 50) / 100, center[1], center[2]] for i in rises]
                        elif axis == 'y':
                            centers = [[center[0], center[1] + width * (i - 50) / 100, center[2]] for i in rises]
                        elif axis == 'z':
                            centers = [[center[0], center[1], center[2] + width * (i - 50) / 100] for i in rises]
                    else:
                        rises = [-1]
                        centers = [center]

                    for ccount, center in zip(rises, centers):
                        if ccount >= 0:
                            fo = f"{kind}-{axis}-{field}-out{out}-rise{ccount}.png"
                        else:
                            fo = f"{kind}-{axis}-{field}-out{out}.png"
                        fo = args.prefix + fo
                        if (not args.overwrite) and os.path.exists(os.path.join(args.to, fo)):
                            print(f"{args.to}/{fo} exists. Skipping")
                            continue
                        if axis in ['face', 'edge']:
                            los = dires[axis][0]
                            north = dires[axis][1]
                        print("Plotting", fo)
                        if kind == 'slc':
                            if axis in ['x', 'y', 'z']:
                                # if args.no_ytfast:
                                # do not use ytfast.SlicePlot
                                if 1:
                                    p = yt.SlicePlot(ds, axis, thefield, center=center, width=width)
                                else:
                                    p = ytfast.SlicePlot(ds, axis, thefield, center=center, width=width,
                                                         force_redo=args.overwrite_cache)
                            else:
                                p = yt.OffAxisSlicePlot(
                                    ds, los, thefield, center=center,
                                    width=width, north_vector=north)
                        elif kind == 'prj':
                            if axis in ['x', 'y', 'z']:
                                if args.no_ytfast:
                                    p = yt.ProjectionPlot(ds, axis, thefield, center=center, width=width,
                                                          weight_field=('gas', 'density'),
                                                          max_level=18,
                                                          )
                                else:
                                    try:
                                        p = ytfast.ProjectionPlot(ds, axis, thefield, center=center, width=width,
                                                                  weight_field=('gas', 'density'),
                                                                  level=18,
                                                                  force_redo=args.overwrite_cache)
                                    except yt.utilities.exceptions.YTPixelizeError:
                                        p = yt.ProjectionPlot(ds, axis, thefield, center=center, width=width,
                                                              weight_field=('gas', 'density'))
                            else:
                                p = yt.OffAxisProjectionPlot(
                                    ds, los, thefield, center=center, width=width,
                                    weight_field=('gas', 'density'), #max_level=l_max,
                                    north_vector=north)
                        elif kind == 'colden':
                            if axis in ['x', 'y', 'z']:
                                use_yt_fast = not args.no_ytfast
                                if use_yt_fast:
                                    try:
                                        if args.den_limit_lower is None:
                                            p = ytfast.ProjectionPlot(ds, axis, thefield, center=center, width=width, force_redo=args.overwrite_cache)
                                        else:
                                            p = ytfast.ProjectionPlot(ds, axis, thefield, center=center, width=width, force_redo=args.overwrite_cache, data_source=high_density_region)
                                    except yt.utilities.exceptions.YTPixelizeError:
                                        use_yt_fast = False # fall back
                                if not use_yt_fast:
                                    if args.den_limit_lower is None:
                                        p = yt.ProjectionPlot(ds, axis, thefield, center=center, width=width)
                                    else:
                                        p = yt.ProjectionPlot(ds, axis, thefield, center=center, width=width, data_source=high_density_region)
                            else:
                                p = yt.OffAxisProjectionPlot(ds, los, thefield, center=center, width=width, north_vector=north)
                        if field in ['den', 'density']:
                            if kind in ['slc', 'prj']:
                                plotutils.den_setup(p)
                            elif kind in ['colden']:
                                plotutils.colden_setup(p)
                            #p.set_colorbar_label(thefield, r"log n [cm$^{-3}$]")
                        # if 'zlim' in args.params:
                        #     if field in args.params['zlim']:
                        #         p.set_zlim(thefield, *args.params['zlim'][field])
                        if thefield[1] == "density" and args.den_zlim is not None:
                            p.set_zlim(thefield, *args.den_zlim)
                        elif thefield[1] == "temperature" and args.T_zlim is not None:
                            p.set_zlim(thefield, *args.T_zlim)
                        elif thefield[1] == "xHII" and args.xHII_zlim is not None:
                            p.set_zlim(thefield, *args.xHII_zlim)
                        elif thefield[1] == "grid_level":
                            # p.set_cmap(thefield, "tab10")
                            p.set_cmap(thefield, "tab20")
                            p.set_log(thefield, False)
                            # p.set_zlim(thefield, 3.5, 13.5)
                            p.set_zlim(thefield, -0.5, 19.5)
                        # overplot sink particles
                        if args.sink:
                            is_id = False
                            lims = [0.1, 10]
                            if args.mass_lim is not None:
                                lims = args.mass_lim
                            if args.sink_color is None:
                                colors = cm.Greens
                            else:
                                colors = cm.get_cmap(args.sink_color)
                            r.overplot_sink_with_id(
                                p, out, center, width/2, is_id=is_id, colors=colors,
                                zorder='mass', withedge=1, lims=lims)
                        # add a rule bar
                        # if args.annotate_scale:
                        #     # p.annotate_scale(max_frac=0.3, min_frac=0.08)
                        #     if args.scale is not None:
                        #         coeff, unit = float(args.scale.split(',')[0]), args.scale.split(',')[1]
                        #         p.annotate_scale(coeff=coeff, unit=unit)
                        #     else:
                        #         p.annotate_scale(max_frac=0.4, min_frac=0.15, )
                        p.set_figure_size(args.figsize)
                        if args.t0 is not None:
                            if args.time_unit is None:
                                thefmt = ".1f"
                            elif args.time_unit == "Myr":
                                thefmt = ".1f"
                            elif args.time_unit == "kyr":
                                thefmt = ".0f"
                            p.annotate_timestamp(time_format='t = {time:' + thefmt + '} {units}',
                                                 time_offset=timeoffset, corner="upper_left",
                                                 time_unit=args.time_unit)

                        if args.contour is not None:
                            p.annotate_contour(thefield, ncont=1, clim=(args.contour, args.contour), label=True, take_log=False, 
                                               plot_args={'colors': 'white', 'linewidths': 1.5})

                        if not args.show_axes:
                            # p.axes.get_xaxis().set_visible(False)
                            # p.axes.get_yaxis().set_visible(False)
                            p.hide_axes(draw_frame=True)
                            p.set_axes_unit('pc')
                        if args.draw_box is not None:
                            plotutils.annotate_box(p, boxcenter, box_width, axis)
                        f = p.export_to_mpl_figure((1,1))
                        ax = f.axes[0]
                        if field in ['den', 'density']:
                            cb = f.axes[1]
                            yticks = cb.get_yticks()
                            lbs = [str(int(log10(x))) for x in yticks]
                            # cb.set_yticklabels(['3', '4', '5', '6', '7', '8', '9'])
                            cb.set_yticklabels(lbs)
                            if kind in ['slc', 'prj']:
                                cb.set_ylabel(r"$\log$ density (cm$^{-3}$)")
                            elif kind in ['colden']:
                                cb.set_ylabel(r"$\log \ \Sigma$ (g cm$^{-2}$)")
                        if (not args.show_axes) and args.annotate_scale:
                            # p.annotate_scale(max_frac=0.3, min_frac=0.08)
                            if args.scale is not None:
                                coeff_str, unit = args.scale.split(',')
                                coeff = float(coeff_str)
                                leng = ramses.to_boxlen((coeff, unit), ds) / ramses.to_boxlen((1, 'pc'), ds)
                                label = f"{coeff_str} {unit}"
                                if args.scaleloc == "lower_left":
                                    left, right = 0.06, None
                                else:
                                    left, right = None, 0.94
                                plotutils.add_scalebar(
                                    ax, leng, label=label, left=left, right=right, 
                                    fontsize="large")
                                # p.annotate_scale(coeff=coeff, unit=unit)
                            else:
                                pass
                                # p.annotate_scale(max_frac=0.4, min_frac=0.15, )
                        # add axis labels
                        if not args.show_axes:
                            plotutils.annotate_axis_label(ax, axis)
                        if args.dry_run:
                            f.savefig(fo, dpi=300)
                            return
                        f.savefig(os.path.join(args.to, fo), dpi=300)
                        print(f"{args.to}/{fo} saved")
                        plt.close('all')
                        del p
                        # p.save((args.to, fo), mpl_kwargs={"dpi": 300})
    return


