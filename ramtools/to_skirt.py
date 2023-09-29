import os
import sys
import numpy as np
from numpy.linalg import norm
import astropy.units as U
import argparse

# from pkgs.ramses import Ramses, NoSinkParticle
# from pkgs.utils import utilities, units
from ramtools.ramses import Ramses, NoSinkParticle
from ramtools import utilities, units


def to_skirt(jobdir, output, l_max, fn_out, family='BC', center=None, width_box_frac=1.0, Z=0.02, letdie=False, skip_exist=True):
    """ Make particle data for SKIRT run

    Args:
        jobdir (str): RAMSES job directory
        output (int): the output number
        l_max (int): maximum level of AMR refinement. This is used to determine the sink particle size, assume R_sink = 2*dx and dx = box_size / 2**l_max.
        fn_out (str): output filename
        family (str): SED family. Default: 'BC'
        center (float or array): center of the box in boxlen units. Default: (.5, .5, .5)
        width_box_frac: width of the box as a fraction of the full simulation box. Default: 1.0
        Z (float): metallicity. Default: solar (0.02)
        letdie (bool): whether to let stars die depending on their ages. Default: False
    """

    if skip_exist and os.path.isfile(fn_out):
        print(f"{fn_out} exists. Skpipped")
        return
    header = "# Star particle for testing SED family\n"\
        "# Column 1: position x (pc)\n"\
        "# Column 2: position y (pc)\n"\
        "# Column 3: position z (pc)\n"\
        "# Column 4: size h (pc)\n"
    if family == 'BB':
        header += "# Column 5: radius (km)\n"
        header += "# Column 6: temperature (K)"
    elif family == 'BC':
        header += "# Column 5: mass (Msun)\n"
        header += "# Column 6: metallicity (1)\n"
        header += "# Column 7: age (Myr)"
    elif family == 'CK':  # Kurucz SED
        header += "# Column 5: radius (km)\n"
        header += "# Column 6: metallicity (1)\n"
        header += "# Column 7: Teff (K)\n"
        header += "# Column 8: g (m/s2)"
    else:
        print('Unrecognized SED', family)
        return 1

    r = Ramses(jobdir=jobdir)

    unitl2pc = r.unit_l_code / units.pc
    boxlen = r.boxlen                   # in unit_l_code
    dx = boxlen / 2**l_max
    if center is None:
        center = 0.5 * boxlen           # scalar
    else:
        center = np.array(center) * boxlen  # array
    output = int(output)
    try:
        loc = r.get_sink_positions(output)  # code unit
    except NoSinkParticle:
        np.savetxt(fn_out, np.array([]), header=header, comments='')
        print("No sink particle file found.", fn_out, 'saved with empty content.')
        return
    loc -= center   # code unit
    width = width_box_frac * boxlen   # code unit
    pick = np.max(np.abs(loc), axis=1) <= width / 2
    if letdie:
        is_alive = r.is_sink_alive(output, mass_shift=0.4)
        pick = np.logical_and(pick, is_alive)
    n = np.sum(pick)
    if n == 0:
        with open(fn_out, 'w') as f:
            f.write(header)
        return
    data = np.zeros([n, 8])
    data[:, :3] = loc[pick, :] * unitl2pc
    # assuming sink particle size is 2*dx
    data[:, 3] = 2 * dx
    mass_to_temp_v = np.vectorize(utilities.mass_to_temp)
    mass_to_radius_v = np.vectorize(utilities.mass_to_radius)
    if family == 'BB':
        # BB
        ncol = 6
        mass = r.get_sink_masses(output)
        T19 = mass_to_temp_v(mass)
        R19 = mass_to_radius_v(mass)
        Rsun = 695700.  # in km
        data[:, 4] = R19[pick] * Rsun
        data[:, 5] = T19[pick]
    elif family == 'BC':
        # BC
        ncol = 7
        data[:, 4] = r.get_sink_masses(output)[pick]
        data[:, 5] = Z
        particles = np.loadtxt(r.get_sink_path(output), delimiter=',')
        data[:, 6] = particles[pick, 11] * r.unit_t * U.s.to('Myr')
    elif family == 'CK':
        ncol = 8
        mass = r.get_sink_masses(output)    # in Msun
        T19 = mass_to_temp_v(mass)
        R19 = mass_to_radius_v(mass)  # in Rsun
        Rsun = 695700.  # in km
        data[:, 4] = R19[pick] * Rsun
        data[:, 5] = Z
        data[:, 6] = T19[pick]
        Rsun2m = 695700000.0  # SI
        G = 6.6743e-11  # SI
        Msun2kg = 1.988409870698051e+30  # SI
        data[:, 7] = G * (mass[pick] * Msun2kg) / (R19[pick] * Rsun2m)**2
    else:
        print('Unrecognized SED', family)
        return 1
    np.savetxt(fn_out, data[:, :ncol], header=header, comments='')
    print(fn_out, 'saved')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('jobid', )
    parser.add_argument('out', type=int, )
    parser.add_argument('of', help='filename of the output particle file')
    parser.add_argument('-c', nargs='+', type=float, default=None,
                        help="""a list representing the center of the particles
                        in code unit (~pc). e.g. -c 0.5 0.5 0.5""")
    parser.add_argument('-s', default='BC', help='type of SED: BB, BC')
    parser.add_argument('-w', default=np.inf, type=float,
                        help='width of the box (in code unit ~pc) inside'\
                        'which the particles are allowed to locate')
    parser.add_argument('-width_boxlen', default=np.inf, type=float,
                        help='width of the box (in box unit) inside'\
                        'which the particles are allowed to locate')
    parser.add_argument('-die', action="store_true", default=False,
                        help='let stars die')

    args = parser.parse_args()

    to_skirt(args.jobid, args.out, args.of, args.s, args.c, args.w,
             args.width_boxlen, letdie=args.die)
