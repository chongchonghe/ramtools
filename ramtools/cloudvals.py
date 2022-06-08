"""
Cloudvals.py generates a text fi with important cloud characteristics at each time step. The energy option calculates thermal, turbulent, magnetic, and gravitational potential energies. The density option calculates cloud mass, volume, and bulk density.

"""

import os
import numpy as np
import yt
import matplotlib.pyplot as plt
import argparse
from yt.units import erg
from ramses import Ramses
from plotutils import quick_plot_prj, den_setup

# yt.enable_parallelism()

OUT = './images/cloudvals'

def main():
    parser = argparse.ArgumentParser(description = 'Plot as a function of time')
    parser.add_argument('jobdir', type = str, nargs='+', help='Directory containing input data')
    parser.add_argument('--frames', type=int, nargs=2, default=[1,10], help='Range of Frames to plot. Defaults to 1-10')
    # parser.add_argument('--out', type = str, nargs = 1, default='.', help='Output location')
    parser.add_argument('--field', type = str, nargs = 1, default = ['energy'], help='Values to calculate, either "energy for energy types, or "density for physical cloud parameters"')
    parser.add_argument('--nocut', action = 'store_true', help='Do not cut on cloud')
    parser.add_argument('--coreonly', action = 'store_true', help='Do not cut on cloud')
    parser.add_argument('--densityandT', action = 'store_true', help='Cut on density (>2 cm-3) and temperature (<1000 K)')
    parser.add_argument('--ref', type=str, default='center',
                        help='choose which pixel to use as the reference to calculate the gravitatoinal potential')
    args = parser.parse_args()
    jobdirs = args.jobdir
    field = args.field[0]
    frames = args.frames
    # out = args.out
    out = OUT
    os.makedirs(out, exist_ok=True)
    nocut = args.nocut
    coreonly = args.coreonly
    densityandT = args.densityandT
    ref = args.ref

    colors = ['tab:blue','tab:green','tab:purple','tab:red']

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()

    for n, jobdir in enumerate(jobdirs):
        times = []
        E_B_list = []
        E_turb_list = []
        E_therm_list = []
        E_pot_list = []
        ram = Ramses(jobdir)
        if nocut:
            cutlabel = '_nocut'
        elif coreonly:
            cutlabel = '_coreonly'
        elif densityandT:
            cutlabel = '_DensityAndT'
        else: # default
            cutlabel = ''

        ref_label = '_' + ref

        jobname = os.path.basename(jobdir)
        doit = False
        if field == 'energy':
            fo = f"{out}/Energy_{jobname}{cutlabel}{ref_label}.txt"
            # if not os.path.exists(fo):
            if 1:
                doit = True
                fi = open(fo, 'w')
                fi.write(f"Time, Magnetic_Energy, Thermal_Energy, Gravitational_Potential, Turbulent_Energy, virial_ratio\n")
                fi.write(f"Myrs, ergs, ergs, ergs, ergs, 1\n")
        elif field == 'density':
            fi = open(f"{out}/Density_{jobname}{cutlabel}.txt", 'a')
            fi.write(f"Time, Cloud_Mass, Cloud_Volume, Average_Density\n")
            fi.write(f"Myrs, cm^3, g, cm^-3\n")
        else:
            print(f"field = {field}")
            return

        if doit:
            for i in range(frames[0],frames[1]):
                try:
                    ds = ram.load_ds(i)
                    time_myrs = ds.current_time.in_units('Myr') #Converts from code units to Megayears
                    times.append(time_myrs)
                except:
                    continue

                ds = ram.load_ds(i)
                mag_unit = ds.magnetic_unit
                alldat = ds.all_data()
        #        mask = reg['density'].in_units("cm**-3", equivalence="number_density", mu=1.4) > 100 #Unused Density Cut
        #        reg = alldat.cut_region(["obj['density'].in_units('cm**-3', equivalence='number_density', mu=1.4) > 100"])
                if nocut:
                    reg = alldat
                elif coreonly:
                    reg = alldat.cut_region(["(obj['temperature'] < 1000) & (obj['density'].in_units('cm**-3', equivalence='number_density', mu=1.4) > 900)"])
                elif densityandT:
                    reg = alldat.cut_region(["(obj['temperature'] < 1000) & (obj['density'].in_units('cm**-3', equivalence='number_density', mu=1.4) > 2)"])
                else: # default
                    reg = alldat.cut_region(["obj['temperature'] < 1000"])

                if field == 'energy':

                    # -- Chong-Chong He
                    # calculate the 'actual' gravitational potential energy of one pixel
                    ad = alldat
                    x = ad['x']
                    y = ad['y']
                    z = ad['z']
                    mass = ad['cell_mass']
                    assert len(x) == len(y) == len(z) == len(mass), \
                        f"len x y z mass = {len(x), len(y), len(z), len(mass)}, shape of particle_pos: {r.shape}"
                    den = ad['density']
                    if ref == 'max':
                        pick = np.argmax(den)
                    elif ref == 'listcenter':
                        pick = len(mass) // 2
                    elif ref == 'center':
                        bx = np.abs(x/ds.length_unit - .5) < 1e-6
                        by = np.abs(y/ds.length_unit - .5) < 1e-6
                        bz = np.abs(z/ds.length_unit - .5) < 1e-6
                        picks = np.logical_and(np.logical_and(bx, by), bz)
                        if np.sum(picks) == 0:
                            print(f"fails to pick the center in frame {i}")
                            continue
                        pick = np.argmax(picks)
                    elif ref == 'corner':
                        loc = [.8, .8, .8]
                        bx = np.abs(x/ds.length_unit - loc[0]) < 1e-2
                        by = np.abs(y/ds.length_unit - loc[1]) < 1e-2
                        bz = np.abs(z/ds.length_unit - loc[2]) < 1e-2
                        picks = np.logical_and(np.logical_and(bx, by), bz)
                        if np.sum(picks) == 0:
                            print(f"fails to pick the center in frame {i}")
                            continue
                        pick = np.argmax(picks)
                    elif ref == 'corner2':
                        loc = [.3, .3, .7]
                        bx = np.abs(x/ds.length_unit - loc[0]) < 1e-2
                        by = np.abs(y/ds.length_unit - loc[1]) < 1e-2
                        bz = np.abs(z/ds.length_unit - loc[2]) < 1e-2
                        picks = np.logical_and(np.logical_and(bx, by), bz)
                        if np.sum(picks) == 0:
                            print(f"fails to pick the center in frame {i}")
                            continue
                        pick = np.argmax(picks)
                    # dist = np.linalg.norm(x - x[pick], axis=1)
                    dist = np.sqrt((x - x[pick])**2 + (y - y[pick])**2 + (z - z[pick])**2)
                    ipot = np.sum(mass[:pick] / dist[:pick]) + np.sum(mass[pick+1:] / dist[pick+1:])
                    ipot *= -1. * yt.units.G  # cgs
                    pot_unit = (ds.velocity_unit)**2 #Unit for potential in cgs
                    pot_ground = ipot - ad['potential'][pick] * pot_unit

                    bx_avg = ((reg['x-Bfield-left']+reg['x-Bfield-right'])/2 * mag_unit)
                    by_avg = ((reg['y-Bfield-left']+reg['y-Bfield-right'])/2 * mag_unit)
                    bz_avg = ((reg['z-Bfield-left']+reg['z-Bfield-right'])/2 * mag_unit)

                    E_B = ((bx_avg**2+by_avg**2+bz_avg**2)/(8*np.pi)*reg['cell_volume']).in_cgs()

                    vx = reg['x-velocity']
                    vy = reg['y-velocity']
                    vz = reg['z-velocity']
                    E_turb = (0.5*((vx*vx)+(vy*vy)+(vz*vz))*reg['cell_mass']).in_cgs()
            #        print(f'Vel = {np.sqrt(np.sum(E_turb)/np.sum(reg["cell_mass"]))}')

        #            pot_unit = (ds.length_unit**2/(ds.time_unit**2)) #Unit for potential in cgs
                    E_pot = ((reg['potential'] * pot_unit + pot_ground) * reg['cell_mass']).in_cgs()

                    press = reg['Pressure']
                    E_therm = (3/2*press*reg['cell_volume']).in_cgs()

                    E_B_sum = np.sum(E_B)
                    E_turb_sum = np.sum(E_turb)
                    E_pot_sum = np.sum(E_pot)
                    E_therm_sum = np.sum(E_therm)
                    E_B_list.append(E_B_sum)
                    E_turb_list.append(E_turb_sum)
                    E_pot_list.append(E_pot_sum)
                    E_therm_list.append(E_therm_sum)
                    alpha = - E_turb_sum / E_pot_sum
        #            print(E_pot_list)
                    fi.write(f"{time_myrs.value}, {E_B_sum.value}, {E_therm_sum.value}, {E_pot_sum.value}, {E_turb_sum.value}, {alpha} \n")
                    print(f"frame = {i}, time = {time_myrs}, E_B = {E_B_sum}, E_t = {E_turb_sum}, E_g = {E_pot_sum}, E_th = {E_therm_sum}")
                    print(f"E_turb / E_pot = {alpha}")

                elif field == 'density':
                    vol = np.sum(reg['cell_volume'])
                    mass = np.sum(reg['cell_mass'])
                    rhoavg = (mass/vol).in_units('cm**-3', equivalence='number_density', mu=1.4)
                    fi.write(f"{time_myrs.value}, {mass.value}, {vol.value}, {rhoavg.value}\n")
            fi.close()

        # make a figure
        with open(fo) as f:
            vals = f.readline().split(',')
            units = f.readline().split(',')
        data = np.loadtxt(fo, skiprows=2, delimiter=',')
        assert len(data) > 0, f'{fo} is empty'
        f, ax = plt.subplots()
        t = data[:, 0]
        for i in range(data.shape[1]-1):
            val = vals[i+1].replace('\n', '')
            sign = 1.
            un = units[i+1].replace('\n', '')
            lb = f"{val} [{un}]"
            if 'Grav' in val:
                ax.plot(t, -1. * data[:, i+1], label='- ' + lb)
                ax.plot(t, -1. * 0.5 * data[:, i+1], label='- 0.5 ' + lb)
            else:
                ax.plot(t, data[:, i+1], label=lb)
        # mark tRelax
        ax.axvline(ram.tRelax, color='gray', ls='--')
        # mark sink formation
        for i in range(frames[0],frames[1]):
            try:
                ram.get_sink_masses(i)
                tsink = ram.get_time(i)
                break
            except:
                pass
        ax.axvline(tsink, color='red', ls='--')
        ax.legend()
        ax.set_ylim([0, None])
        plt.savefig(f"{out}/fig_energy_{jobname}{cutlabel}{ref_label}.png", dpi=300)

        f, ax = plt.subplots()
        t = data[:, 0]
        for i in range(data.shape[1]-1):
            val = vals[i+1].replace('\n', '')
            if 'Grav' in val:
                grav = data[:, i+1]
                continue
            if 'Turb' in val:
                turb = data[:, i+1]
                continue
        ax.plot(t, -1. * turb / grav, label='e_K / |e_G|' )
        ax.plot(t, -1. * 2 * turb / grav, label='e_K / 0.5 |e_G|' )
        ax.axhline(0.5, color='k', ls='--')
        # mark tRelax
        ax.axvline(ram.tRelax, color='gray', ls='--')
        # mark sink formation
        for i in range(frames[0],frames[1]):
            try:
                ram.get_sink_masses(i)
                tsink = ram.get_time(i)
                break
            except:
                pass
        ax.axvline(tsink, color='red', ls='--')
        ax.legend()
        ax.set_ylim([0, None])
        plt.savefig(f"{out}/fig_virial_{jobname}{cutlabel}{ref_label}.png", dpi=300)

