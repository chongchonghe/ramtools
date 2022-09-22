"""
Calculate some parameters of a cloud at given snapshot(s).
"""
import os
import numpy as np
import yt
import yt.units as unyt
import json
import ramtools as rt
from . import utilities as util
from . import plotutils as pu
from .units import mH
# from .cacherun import CacheRun, RamsesSnapshot
from .cacherun import CacheRun


def get_cloud_value(ds, var, dataset=None, filter=None, pot_ref='center',
                    pot_test_width=2**-12):
    """
    Return the sum of a given quantity (e.g. the total mass.) over a region defined by filter of a given cloud (identified by ds). 

    Args:
        ds: yt.ds, instance of ds.load(...)
        var: the variable name, one of totmass, totvolume, kin_ene,
            th_ene, mag_ene, pot_ene.
        dataset: (optional) the dataset, e.g. ds.sphere(...)
        filter: (string) filter to the dataset
        pot_ref: the reference point of calculating potential energy. Refer
            to the code
        pot_test_width: Refer to the code.

    Returns:
        YT float.
    """

    if dataset is None:
        # ad = ds.all_data()
        ad = ds.sphere('c', radius=0.5)
    else:
        ad = dataset
    if filter is None:
        reg = ad
    else:
        reg = ad.cut_region([filter])
    if var == "totmass":
        return reg.quantities.total_quantity(
            [("gas", "cell_mass")]
        )
    if var == "totvolume":
        return reg.quantities.total_quantity(
            [("gas", "cell_volume")]
        )
    if var == 'kin_ene':
        vx = reg['x-velocity']
        vy = reg['y-velocity']
        vz = reg['z-velocity']
        E_turb = (0.5 * ((vx * vx) + (vy * vy) + (vz * vz)) * reg[
            'cell_mass']).sum()
        return E_turb
    if var == 'therm_ene':
        press = reg['Pressure']
        E_therm = (3 / 2 * press * reg['cell_volume']).sum()
        return E_therm
    if var == 'mag_ene':
        # unit_B = util.get_unit_B(ds)
        unit_B = util.get_unit_B_new(ds)    # CORRECTION on 2022-9-15 by CCH
        bx_avg = (reg['x-Bfield-left'] + reg[
            'x-Bfield-right']) / 2 * unit_B
        by_avg = (reg['y-Bfield-left'] + reg[
            'y-Bfield-right']) / 2 * unit_B
        bz_avg = (reg['z-Bfield-left'] + reg[
            'z-Bfield-right']) / 2 * unit_B
        E_B = ((bx_avg ** 2 + by_avg ** 2 + bz_avg ** 2) / (8 * np.pi) *
               reg['cell_volume']).sum()
        return E_B
    # end if var == 'mag_ene':
    if var == 'pot_ene':
        dx = ds.length_unit / 2 ** 12  # the size of the box to exclude
        dx_cm = float(dx)
        lu_cm = float(ds.length_unit)
        if isinstance(pot_ref, str):
            if pot_ref == 'max':
                # pick = np.argmax(den)
                # bx = np.abs(x / ds.length_unit - .5) < 1e-5
                # by = np.abs(y / ds.length_unit - .5) < 1e-5
                # bz = np.abs(z / ds.length_unit - .5) < 1e-5
                # picks = np.logical_and(np.logical_and(bx, by), bz)
                # if np.sum(picks) == 0:
                #     print(f"fails to pick the center in frame {i}")
                # pick = np.argmax(picks)
                xt, ytt, zt, mt, pott = ad.argmax(
                    'density', axis=['x', 'y', 'z', 'cell_mass', 'potential'])
            elif pot_ref == 'bottom_left':
                p1 = ad.include_inside('x', 0, lu_cm * pot_test_width)
                p2 = p1.include_inside('y', 0, lu_cm * pot_test_width)
                p3 = p2.include_inside('z', 0, lu_cm * pot_test_width)
                xt = p3['x'][0]
                ytt = p3['y'][0]
                zt = p3['z'][0]
                # mt = p3['cell_mass'][0]
                pott = p3['potential'][0]
            elif pot_ref == 'bottom_right':
                p1 = ad.include_inside('x', lu_cm - lu_cm * pot_test_width,
                                       lu_cm)
                p2 = p1.include_inside('y', 0, lu_cm * pot_test_width)
                p3 = p2.include_inside('z', 0, lu_cm * pot_test_width)
                xt = p3['x'][0]
                ytt = p3['y'][0]
                zt = p3['z'][0]
                # mt = p3['cell_mass'][0]
                pott = p3['potential'][0]
            else:
                print("Unrecognized pot_ref:", pot_ref)
                values.append(np.nan)
        else:
            xt, ytt, zt = np.array(pot_ref) * ds.length_unit
            point_obj = ds.point(np.array(pot_ref) * ds.length_unit)
            pott = point_obj['potential'][0]
        ep1 = ad.exclude_inside('x', float(xt - dx), float(xt + dx))
        ep2 = ad.exclude_inside('y', float(ytt - dx), float(ytt + dx))
        ep3 = ad.exclude_inside('z', float(zt - dx), float(zt + dx))
        x = ep3['x']
        y = ep3['y']
        z = ep3['z']
        m = ep3['cell_mass']
        dist = np.sqrt(
            (x - xt) ** 2 + (y - ytt) ** 2 + (z - zt) ** 2)
        ipot = -1. * unyt.gravitational_constant * np.sum(m / dist)
        try:
            pot_unit = 1.
            pot_ground = ipot - pott * pot_unit
        except unyt.exceptions.UnitOperationError:
            pot_unit = (ds.velocity_unit) ** 2  # Unit for potential in cgs
            pot_ground = ipot - pott * pot_unit
        totpot = ((reg['potential'] * pot_unit + pot_ground) * reg[
            'cell_mass']).sum()
        return totpot
    # end if var == 'pot_ene'
    return np.nan


def describe_cloud(jobfolder, out, variables, filter):
    """
    Describe the parameters of a cloud at a given snapshot

    Args:
        jobfolder:
        out:
        variables (list of strings): a list of variables to calculate
        filter: gas region filter

    Returns:
        A dictionary of all parameters
    """

    # snap = RamsesSnapshot(jobfolder, out)
    # ds = snap._ds
    # dataset = ds.sphere('c', radius=0.5)
    # filter = "obj['density'].in_units('cm**-3', equivalence='number_" \
    #          "density', mu=1.4) > 3"
    ds = rt.utilities.my_yt_load(jobfolder, out)
    # nthresh = 3 * 1.4 * mH
    # filter = f"obj['density'] > {nthresh}"
    ret = {}
    for var in variables:
        # y = get_cloud_value(ds, var, dataset)
        # task = CacheRun(get_cloud_value)
        # task.ForceReplaceCache()
        # y = task.Run(snap, var, filter=filter)
        y = CacheRun(get_cloud_value, 1)(ds, var, filter=filter)
        ret[var] = (float(y.in_cgs()), str(y.in_cgs().units))
    # make slice plot of the density and temperature
    return ret


def describe_cloud2(jobfolder, out, variables, slicepath, fo=None):
    """
    Describe the parameters of a cloud at a given snapshot

    Args:
        jobfolder:
        out:
        variables (list of strings): a list of variables to calculate
        fo (str): name of the output json file
        dataset: YT dataset, e.g. ds.all_data()
        filter: cut_region filter
        pot_ref:
        pot_test_width:

    Returns:
        None
    """

    # snap = RamsesSnapshot(jobfolder, out)
    # ds = snap._ds
    # dataset = ds.sphere('c', radius=0.5)
    # filter = "obj['density'].in_units('cm**-3', equivalence='number_" \
    #          "density', mu=1.4) > 3"
    ds = rt.utilities.my_yt_load(jobfolder, out)
    nthresh = 3 * 1.4 * mH
    filter = f"obj['density'] > {nthresh}"
    ret = {}
    for var in variables:
        # y = get_cloud_value(ds, var, dataset)
        # task = CacheRun(get_cloud_value)
        # task.ForceReplaceCache()
        # y = task.Run(snap, var, filter=filter)
        y = CacheRun(get_cloud_value)(ds, var, filter=filter)
        ret[var] = (float(y.in_cgs()), str(y.in_cgs().units))
    # make slice plot of the density and temperature
    os.makedirs(slicepath, exist_ok=True)
    for axis in ['x', 'y', 'z']:
        for field in [('gas', 'density'), ('gas', 'temperature')]:
            figo = f"{slicepath}/{out}_{axis}_{field[1]}.png"
            if os.path.exists(figo):
                continue
            try:
                print(reg)
            except:
                dataset = ds.sphere('c', radius=0.5)
                reg = dataset.cut_region([filter])
            p = yt.SlicePlot(ds, axis, field, data_source=reg)
            if field[1] == 'density':
                pu.den_setup(p)
            elif field[1] == 'temperature':
                pu.T_setup(p)
            p.save(figo)
    if fo is None:
        return ret
    else:
        with open(fo, 'w') as f:
            json.dump(ret, f, indent=2)
        return


def check_pot_ene():

    # rp = rt.Ramses("tests/Job2")
    # for pot_ref in ['bottom_left', 'bottom_right', 'max']:
    #     ret = rt.ramses.get_cloud_values(rp.jobPath, 40, ['kin_ene', 'pot_ene'],
    #                                      pot_ref=pot_ref, pot_test_width=2**-7)
    #     print(ret)
    print("start")
    snap = rt.cacherun.RamsesSnapshot("tests/Job2", 40)
    rets = []
    np.random.seed(2021)
    ran1 = np.random.random(3)
    ran2 = np.random.random(3)
    ran3 = np.random.random(3)
    reflist = ['bottom_left', 'bottom_right', 'max', ran1, ran2, ran3]
    for pot_ref in reflist:
        ret = rt.cacherun.CacheRun(get_cloud_values)(
            snap, ['kin_ene', 'pot_ene'], pot_ref=pot_ref, pot_test_width=2**-7)
        rets.append(ret)
    for i, j in zip(rets, reflist):
        print(j)
        print(i)
    print("Are they close?")

    # Will to work. Can't pickle local object 'PlotWindow'
    # # def proj_ad_wrapper(ds, *args, **kwargs):
    # #     return yt.ProjectionPlot(ds.all_data(), *args, **kwargs)
    # p = hash2.HashRun(yt.ProjectionPlot).Run(
    #     ds, 'x', ('gas', 'density'))
    # p.save('t.png')
    return
