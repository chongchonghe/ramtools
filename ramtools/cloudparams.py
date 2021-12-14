"""
Calculate some parameters of a cloud at given snapshot(s).
"""
import numpy as np
import yt.units as unyt
from . import utilities as util

def get_cloud_values(ds, vars, filter=None, pot_ref='center',
                     pot_test_width=2**-12):

    # ds = yt.load(f"{jobdir}/output_{out:05d}/info_{out:05d}.txt",
    #              fields=util.FIELDS)
    unit_B = util.get_unit_B(ds)
    ad = ds.all_data()
    if filter is None:
        reg = ad
    else:
        reg = ad.cut_region([filter])
    values = []
    for var in vars:
        if var == 'kin_ene':
            vx = reg['x-velocity']
            vy = reg['y-velocity']
            vz = reg['z-velocity']
            E_turb = (0.5 * ((vx * vx) + (vy * vy) + (vz * vz)) * reg[
                'cell_mass']).sum()
            values.append(E_turb.in_cgs())
        elif var == 'therm_ene':
            press = reg['Pressure']
            E_therm = (3/2 * press * reg['cell_volume']).sum()
            values.append(E_therm.in_cgs())
        elif var == 'mag_ene':
            bx_avg = (reg['x-Bfield-left'] + reg[
                'x-Bfield-right']) / 2 * unit_B
            by_avg = (reg['y-Bfield-left'] + reg[
                'y-Bfield-right']) / 2 * unit_B
            bz_avg = (reg['z-Bfield-left'] + reg[
                'z-Bfield-right']) / 2 * unit_B
            E_B = ((bx_avg ** 2 + by_avg ** 2 + bz_avg ** 2) / (8 * np.pi) *
                   reg['cell_volume']).sum()
            values.append(E_B.in_cgs())
        elif var == 'pot_ene':
            dx = ds.length_unit / 2 ** 12   # the size of the box to exclude
            dx_cm = float(dx)
            lu_cm = float(ds.length_unit)
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
                mt = p3['cell_mass'][0]
                pott = p3['potential'][0]
            elif pot_ref == 'bottom_right':
                p1 = ad.include_inside('x', lu_cm - lu_cm * pot_test_width,
                                       lu_cm)
                p2 = p1.include_inside('y', 0, lu_cm * pot_test_width)
                p3 = p2.include_inside('z', 0, lu_cm * pot_test_width)
                xt = p3['x'][0]
                ytt = p3['y'][0]
                zt = p3['z'][0]
                mt = p3['cell_mass'][0]
                pott = p3['potential'][0]
            else:
                print("Unrecognized pot_ref:", pot_ref)
                values.append(np.nan)
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
            values.append(totpot.in_cgs())
        else:
            values.append(np.nan)

    return values


