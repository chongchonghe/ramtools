"""
Cloudvals.py generates a text file with important cloud characteristics at each time step. The energy option calculates thermal, turbulent, magnetic, and gravitational potential energies. The density option calculates cloud mass, volume, and bulk density.

"""

import numpy as np
import yt
import matplotlib.pyplot as plt
import argparse
from ramtools import ramses
from ramtools.plotutils import quick_plot_prj, den_setup
from yt.units import erg

parser = argparse.ArgumentParser(description = 'Plot as a function of time')
parser.add_argument('jobdir', type = str, nargs='+', help='Directory containing input data')
parser.add_argument('--frames', type = int, nargs =2, default= [1,10], help='Range of Frames to plot. Defaults to 1-10')
parser.add_argument('--out', type = str, nargs = 1, default='.', help='Output location')
parser.add_argument('--field', type = str, nargs = 1, default = 'energy', help='Values to calculate, either "energy for energy types, or "density for physical cloud parameters"')
parser.add_argument('--nocut', action = 'store_true', help='Do not cut on cloud')
parser.add_argument('--coreonly', action = 'store_true', help='Do not cut on cloud')
args = parser.parse_args()
jobdirs = args.jobdir
field = args.field[0]
frames = args.frames
out = args.out
nocut = args.nocut
coreonly = args.coreonly

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
    ram = ramses.Ramses(jobdir)
    if nocut:
        cutlabel = '_nocut'
    elif coreonly:
        cutlabel = '_coreonly'
    else:
        cutlabel = ''

    if field == 'energy':
        file = open(f"{out}/Energy_{jobdir[-10:-1]}{cutlabel}.txt", 'a')
        file.write(f"Time, Magnetic_Energy, Thermal_Energy, Gravitational_Potential, Turbulent_Energy\n")
        file.write(f"Myrs, ergs, ergs, ergs, ergs\n")
    elif field == 'density':
        file = open(f"{out}/Density_{jobdir[-10:-1]}{cutlabel}.txt", 'a')
        file.write(f"Time, Cloud_Mass, Cloud_Volume, Average_Density\n")
        file.write(f"Myrs, cm^3, g, cm^-3\n")
        
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
        else:
            reg = alldat.cut_region(["obj['temperature'] < 1000"])
        
        if field == 'energy':
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
            pot_unit = (ds.velocity_unit)**2 #Unit for potential in cgs
            E_pot = (reg['potential'] * pot_unit * reg['cell_mass']).in_cgs()

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
#            print(E_pot_list)
            file.write(f"{time_myrs.value}, {E_B_sum.value}, {E_therm_sum.value}, {E_pot_sum.value}, {E_turb_sum.value} \n")
            
            
            
        
            print(f"E_B = {E_B_list}, E_t = {E_turb_list}, E_g = {E_pot_list}, E_th = {E_therm_list}")
            
        elif field == 'density':
            vol = np.sum(reg['cell_volume'])
            mass = np.sum(reg['cell_mass'])
            rhoavg = (mass/vol).in_units('cm**-3', equivalence='number_density', mu=1.4)
            file.write(f"{time_myrs.value}, {mass.value}, {vol.value}, {rhoavg.value}\n")
    file.close()
    
    
    
